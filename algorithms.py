from __future__ import annotations

import copy
import csv
import json
import math
import warnings
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, List, Sequence, Tuple, TYPE_CHECKING

# 惰性依赖：仅在需要时再导入第三方库，避免默认运行的额外负担
def _get_dilipred_module():
    try:
        import importlib

        return importlib.import_module("dilipred")
    except Exception:
        return None

try:
    from pymoo.core.problem import ElementwiseProblem as _ElementwiseProblem
except Exception:  # pragma: no cover - 缺少依赖时的占位
    class _ElementwiseProblem:  # type: ignore
        def __init__(self, *args, **kwargs):
            raise RuntimeError("pymoo 未安装，进化算法相关功能不可用。")

if TYPE_CHECKING:  # 仅类型检查用途
    import numpy as np  # type: ignore

CONFIG_PATH = Path(__file__).resolve().parent / "config" / "algorithms.json"
_DEFAULT_ALGO_CONFIG = {
    "ga": {"generations": 25, "population_size": 32, "toxicity_weight": 1.0},
    "nsga2": {"generations": 30, "population_size": 40},
    "pso": {
        "iterations": 40,
        "swarm_size": 24,
        "options": {"c1": 1.2, "c2": 1.2, "w": 0.6},
    },
    "aci": {
        "lambda_weight": 0.3,
        "complement_variance_max": 0.25,
        "coverage_floor": 0.0,
    },
}
_ALGO_CONFIG_CACHE: dict | None = None
_DILIPRED_PREDICTOR: Any | None = None


def _normalize_probability(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return _normalize_probability(value[0])
    columns = getattr(value, "columns", None)
    if columns is not None and "value" in columns:
        try:
            series = value["value"]
            first = series.iloc[0] if hasattr(series, "iloc") else series[0]
            return float(first)
        except Exception:
            pass
    if isinstance(value, dict):
        for key in ("probability", "prob", "score", "risk", "dili", "value"):
            if key in value:
                return float(value[key])
    return None


def _load_dilipred_predictor_instance() -> Any | None:
    global _DILIPRED_PREDICTOR
    if _DILIPRED_PREDICTOR is not None:
        return _DILIPRED_PREDICTOR
    dilipred = _get_dilipred_module()
    if dilipred is None:
        return None
    for attr in ("DILIPRedictor", "DILIPred", "DILIPredictor", "Predictor", "Model"):
        factory = getattr(dilipred, attr, None)
        if callable(factory):
            try:
                predictor = factory()
            except TypeError:
                continue
            _DILIPRED_PREDICTOR = predictor
            return predictor
    return None


def _predict_with_dilipred(smiles: str) -> float | None:
    if not smiles:
        return None

    # 直接使用 dilipred 模块所暴露的 Python API；若失败再退化为启发式。
    dilipred = _get_dilipred_module()
    for attr in ("predict_smiles", "predict_dili", "predict"):
        predict_fn = getattr(dilipred, attr, None) if dilipred else None
        if not callable(predict_fn):
            continue
        try:
            result = predict_fn([smiles])
        except TypeError:
            result = predict_fn(smiles)
        score = _normalize_probability(result)
        if score is not None:
            return score

    predictor = _load_dilipred_predictor_instance()
    if predictor is not None:
        for attr in ("predict_proba", "predict"):
            predict_fn = getattr(predictor, attr, None)
            if not callable(predict_fn):
                continue
            try:
                result = predict_fn([smiles])
            except TypeError:
                result = predict_fn(smiles)
            score = _normalize_probability(result)
            if score is not None:
                return score

    return None


def _heuristic_smiles_risk(smiles: str) -> float:
    hetero = sum(1 for ch in smiles if ch in "NOSPF")
    aromatic = smiles.count("c")
    halogen = sum(1 for token in ("Cl", "Br", "I") if token in smiles)
    ring = smiles.count("1")
    risk = 0.25 + 0.02 * hetero + 0.01 * aromatic + 0.05 * halogen + 0.01 * ring
    return float(max(0.0, min(1.0, risk)))


def estimate_hepatotoxicity_from_smiles(smiles: str) -> float:
    """优先使用 DILIPred 估算肝毒性，若不可用则退化为启发式推断。"""

    if not smiles:
        raise ValueError("SMILES 不能为空；无法估算肝毒性。")
    score = _predict_with_dilipred(smiles)
    if score is not None:
        return float(max(0.0, min(1.0, score)))
    warnings.warn(
        "未能通过 DILIPred 获取肝毒性，改用启发式估计（准确度有限）。",
        RuntimeWarning,
        stacklevel=2,
    )
    return _heuristic_smiles_risk(smiles)


def _deep_update(base: dict, overrides: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_update(result[key], value)
        else:
            result[key] = value
    return result


def _load_algorithm_config() -> dict:
    global _ALGO_CONFIG_CACHE
    if _ALGO_CONFIG_CACHE is not None:
        return _ALGO_CONFIG_CACHE
    config = copy.deepcopy(_DEFAULT_ALGO_CONFIG)
    try:
        with CONFIG_PATH.open("r", encoding="utf-8") as fp:
            user_config = json.load(fp)
            config = _deep_update(config, user_config)
    except FileNotFoundError:
        pass
    _ALGO_CONFIG_CACHE = config
    return config


def _get_aci_settings() -> dict:
    """便捷访问 ACI 相关配置。"""

    return _load_algorithm_config().get("aci", _DEFAULT_ALGO_CONFIG["aci"])


@dataclass(frozen=True)
class Ingredient:
    """描述单个中药成分的静态属性。

    Attributes:
        name: 成分名称。
        smiles: 该成分的 SMILES 表达式，用于必要时推断肝毒性。
        hepatotoxicity_score: 外部模型/实验给出的肝毒性打分（越高越差）。
        ob: Oral Bioavailability，若缺失可置 None。
        dl: Drug Likeness，若缺失可置 None。
        mw: 分子量 (MW)。
        alogp: ALogP（脂溶性）。
        h_don: 氢键供体数量。
        h_acc: 氢键受体数量。
        caco2: Caco-2 渗透性（可用 logPapp）。
        bbb: 血脑屏障穿透（logBB）。
        fasa: FASA- 极性表面积比。
        hl: 半衰期（小时）。
    """

    name: str
    smiles: str | None = None
    hepatotoxicity_score: float | None = None
    ob: float | None = None
    dl: float | None = None
    mw: float | None = None
    alogp: float | None = None
    h_don: float | None = None
    h_acc: float | None = None
    caco2: float | None = None
    bbb: float | None = None
    fasa: float | None = None
    hl: float | None = None

    def __post_init__(self) -> None:
        if self.hepatotoxicity_score is not None:
            object.__setattr__(
                self, "hepatotoxicity_score", float(self.hepatotoxicity_score)
            )
            return
        if not self.smiles:
            raise ValueError(
                f"成分 {self.name} 缺少 hepatotoxicity_score，且未提供 SMILES。"
            )
        estimated = estimate_hepatotoxicity_from_smiles(self.smiles)
        object.__setattr__(self, "hepatotoxicity_score", float(estimated))


# ---------------------------
# 数据读取与解析接口
# ---------------------------

_COL_MAP = {
    "name": {"name", "ingredient", "compound", "drug", "名称", "成分"},
    "smiles": {"smiles", "smile"},
    "hepatotoxicity_score": {
        "hepatotoxicity_score",
        "hepatotoxicity",
        "dili",
        "toxicity",
        "toxicity_score",
    },
    "ob": {"ob", "oralbioavailability", "bioavailability"},
    "dl": {"dl", "druglikeness", "drug_likeness"},
    "mw": {"mw", "molecularweight", "molecular_weight"},
    "alogp": {"alogp", "logp", "a_logp"},
    "h_don": {"h_don", "hdon", "hbd", "hbond_donor", "h_donor"},
    "h_acc": {"h_acc", "hacc", "hba", "hbond_acceptor", "h_acceptor"},
    "caco2": {"caco2", "caco-2", "logpapp", "papp"},
    "bbb": {"bbb", "logbb", "bloodbrainbarrier"},
    "fasa": {"fasa", "fasa-"},
    "hl": {"hl", "half_life", "half-life", "t_half"},
}


def _norm_key(key: str) -> str:
    return "".join(ch for ch in key.strip().lower() if ch.isalnum() or ch == "_")


def _map_header_to_attr(header: str) -> str | None:
    nk = _norm_key(header)
    for attr, names in _COL_MAP.items():
        if nk in names:
            return attr
    return None


def _to_float_or_none(v: Any) -> float | None:
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if s == "" or s.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load_ingredients_csv(path: str | Path, encoding: str = "utf-8") -> list[Ingredient]:
    """从 CSV 读取药物成分。

    至少包含 `name`，且 `smiles` 与 `hepatotoxicity_score` 至少其一存在。
    其他可选列：ob, dl, mw, alogp, h_don, h_acc, caco2, bbb, fasa, hl。
    支持常见同义列名（不区分大小写，自动规范化）。
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"未找到文件：{path}")
    with path.open("r", encoding=encoding, newline="") as fp:
        sample = fp.read(4096)
        fp.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except Exception:
            dialect = csv.excel
        reader = csv.DictReader(fp, dialect=dialect)
        if not reader.fieldnames:
            return []
        rows: list[Ingredient] = []
        for raw in reader:
            data: dict[str, Any] = {}
            for key, value in raw.items():
                attr = _map_header_to_attr(key)
                if attr is None:
                    continue
                if attr in {"name", "smiles"}:
                    data[attr] = (value or "").strip() or None
                elif attr == "bbb":
                    fv = _to_float_or_none(value)
                    if fv is None:
                        sv = str(value).strip().lower()
                        if sv in {"y", "yes", "true", "1", "+"}:
                            fv = 1.0
                        elif sv in {"n", "no", "false", "0", "-"}:
                            fv = 0.0
                    data[attr] = fv
                else:
                    data[attr] = _to_float_or_none(value)
            name = (data.get("name") or "").strip()
            if not name:
                continue
            if not data.get("smiles") and data.get("hepatotoxicity_score") is None:
                raise ValueError(
                    f"记录缺少必需信息（必须提供 smiles 或 hepatotoxicity_score）：{name}"
                )
            rows.append(
                Ingredient(
                    name=name,
                    smiles=data.get("smiles"),
                    hepatotoxicity_score=data.get("hepatotoxicity_score"),
                    ob=data.get("ob"),
                    dl=data.get("dl"),
                    mw=data.get("mw"),
                    alogp=data.get("alogp"),
                    h_don=data.get("h_don"),
                    h_acc=data.get("h_acc"),
                    caco2=data.get("caco2"),
                    bbb=data.get("bbb"),
                    fasa=data.get("fasa"),
                    hl=data.get("hl"),
                )
            )
        return rows


def load_ingredients_json(path: str | Path, encoding: str = "utf-8") -> list[Ingredient]:
    """从 JSON 读取药物成分（数组形式，每项为对象）。"""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"未找到文件：{path}")
    data = json.loads(path.read_text(encoding=encoding))
    if not isinstance(data, list):
        raise ValueError("JSON 根应为数组（list）。")
    results: list[Ingredient] = []
    for obj in data:
        if not isinstance(obj, dict):
            continue
        def get_any(keys: Iterable[str]):
            for k in keys:
                if k in obj and obj[k] not in (None, ""):
                    return obj[k]
            return None

        name = get_any(_COL_MAP["name"]) or ""
        smiles = get_any(_COL_MAP["smiles"]) or None
        score = _to_float_or_none(get_any(_COL_MAP["hepatotoxicity_score"]))
        if not name:
            continue
        if not smiles and score is None:
            raise ValueError(
                f"记录缺少必需信息（必须提供 smiles 或 hepatotoxicity_score）：{name}"
            )
        results.append(
            Ingredient(
                name=name,
                smiles=smiles,
                hepatotoxicity_score=score,
                ob=_to_float_or_none(get_any(_COL_MAP["ob"])),
                dl=_to_float_or_none(get_any(_COL_MAP["dl"])),
                mw=_to_float_or_none(get_any(_COL_MAP["mw"])),
                alogp=_to_float_or_none(get_any(_COL_MAP["alogp"])),
                h_don=_to_float_or_none(get_any(_COL_MAP["h_don"])),
                h_acc=_to_float_or_none(get_any(_COL_MAP["h_acc"])),
                caco2=_to_float_or_none(get_any(_COL_MAP["caco2"])),
                bbb=_to_float_or_none(get_any(_COL_MAP["bbb"])),
                fasa=_to_float_or_none(get_any(_COL_MAP["fasa"])),
                hl=_to_float_or_none(get_any(_COL_MAP["hl"])),
            )
        )
    return results


def load_ingredients(path: str | Path, fmt: str | None = None, encoding: str = "utf-8") -> list[Ingredient]:
    """按后缀或指定格式读取药物成分（CSV/JSON）。"""
    path = Path(path)
    ext = path.suffix.lower().lstrip(".")
    format_to_use = (fmt or ("csv" if ext in {"csv", "tsv"} else "json")).lower()
    if format_to_use == "csv":
        return load_ingredients_csv(path, encoding=encoding)
    if format_to_use == "json":
        return load_ingredients_json(path, encoding=encoding)
    raise ValueError(f"不支持的格式：{fmt!r} / 扩展名：.{ext}")


@dataclass(frozen=True)
class FeatureWindow:
    """用于定义 ADME 特征的“宜居区”。

    note 会描述该特征的单位或取值口径，方便事后定位。
    """

    lower: float
    upper: float
    softness: float
    note: str = ""


ADME_FEATURE_ATTRS: dict[str, str] = {
    "mw": "mw",
    "alogp": "alogp",
    "hdon": "h_don",
    "hacc": "h_acc",
    "ob": "ob",
    "caco2": "caco2",
    "bbb": "bbb",
    "dl": "dl",
    "fasa": "fasa",
    "hl": "hl",
}

FEATURE_WINDOWS: dict[str, FeatureWindow] = {
    "mw": FeatureWindow(180.0, 500.0, 120.0, note="Dalton"),
    "alogp": FeatureWindow(-0.5, 5.0, 1.5, note="ALogP"),
    "hdon": FeatureWindow(0.0, 5.0, 2.0, note="count"),
    "hacc": FeatureWindow(0.0, 10.0, 3.0, note="count"),
    "ob": FeatureWindow(30.0, 70.0, 15.0, note="OB%"),
    "caco2": FeatureWindow(0.4, 1.5, 0.3, note="log10(Papp, cm/s)"),
    "bbb": FeatureWindow(-1.0, 1.0, 0.5, note="logBB or mapped分类"),
    "dl": FeatureWindow(0.18, 0.5, 0.1, note="DL"),
    "fasa": FeatureWindow(0.2, 0.6, 0.15, note="FASA-"),
    "hl": FeatureWindow(3.0, 12.0, 3.0, note="hours"),
}

ACI_LAMBDA = 0.3
_SAFE_EPS = 1e-6
_TOTAL_FEATURES = len(ADME_FEATURE_ATTRS)


BITS_PER_BYTE = 8


def _pack_bools(values: Sequence[bool]) -> bytearray:
    data = bytearray(math.ceil(len(values) / BITS_PER_BYTE))
    for idx, value in enumerate(values):
        if value:
            data[idx // BITS_PER_BYTE] |= 1 << (idx % BITS_PER_BYTE)
    return data


def _bit_is_set(bits: bytearray, idx: int) -> bool:
    return bool(bits[idx // BITS_PER_BYTE] & (1 << (idx % BITS_PER_BYTE)))


@dataclass
class CandidateSolution:
    """混合染色体：布尔选择 + 实数配比（使用位存储减少内存占用）。"""

    select_bits: bytearray
    proportions: array
    length: int

    @classmethod
    def from_components(
        cls, selects: Sequence[bool], proportions: Sequence[float]
    ) -> "CandidateSolution":
        return cls(_pack_bools(selects), array("f", proportions), len(selects))

    def iter_selects(self) -> Iterable[bool]:
        for idx in range(self.length):
            yield _bit_is_set(self.select_bits, idx)

    def iter_selected_indices(self) -> Iterable[int]:
        for idx, sel in enumerate(self.iter_selects()):
            if sel:
                yield idx

    def normalized_proportions(self) -> array:
        selects = list(self.iter_selects())
        total = sum(prop for sel, prop in zip(selects, self.proportions) if sel)
        normalized = array("f", [0.0] * self.length)
        if total <= 0:
            selected_indices = [idx for idx, sel in enumerate(selects) if sel]
            if not selected_indices:
                return normalized
            default = 1.0 / len(selected_indices)
            for idx in selected_indices:
                normalized[idx] = default
            return normalized
        for idx, (sel, prop) in enumerate(zip(selects, self.proportions)):
            if sel:
                normalized[idx] = prop / total
        return normalized

    def with_normalized(self) -> "CandidateSolution":
        return CandidateSolution(
            bytearray(self.select_bits), self.normalized_proportions(), self.length
        )


def ensure_selection(
    selects: List[bool], ingredients: Sequence[Ingredient]
) -> List[bool]:
    """保证至少选中一个成分；若全为 False，则按输入顺序启用首个成分以保持确定性。"""
    if any(selects) or not ingredients:
        return list(selects)
    copy = list(selects)
    copy[0] = True
    return copy


def vector_to_candidate(
    vector: Sequence[float], ingredients: Sequence[Ingredient]
) -> CandidateSolution:
    """将 pymoo/pyswarms 的 [0,1] 向量映射为混合染色体。

    Args:
        vector: 兼具“是否选中”和“配比”信息的连续向量。
        ingredients: 组成向量所对应的成分列表，长度应为 vector 长度的一半。

    Returns:
        归一化后的 CandidateSolution。
    """
    # 向量前半段存储“是否选择”布尔值，后半段为连续配比。
    half = len(vector) // 2
    selects = [float(value) > 0.5 for value in vector[:half]]
    selects = ensure_selection(selects, ingredients)
    proportions = [float(value) for value in vector[half:]]
    # 限制配比在 [0.01, 1.0]，防止极端值导致归一化数值漂移。
    proportions = [max(0.01, min(1.0, value)) for value in proportions]
    return CandidateSolution.from_components(selects, proportions)


def _window_score(value: float, window: FeatureWindow) -> float:
    """根据 FeatureWindow 计算单一特征的得分。"""

    if window.lower <= value <= window.upper:
        return 1.0
    delta = window.lower - value if value < window.lower else value - window.upper
    return math.exp(-((delta / max(window.softness, _SAFE_EPS)) ** 2))


def _prepare_feature_value(feature: str, value: float | None) -> float | None:
    """根据特征类型调整原始值（如 BBB 的 0/1 映射）。"""

    if value is None:
        return None
    if feature == "bbb":
        if isinstance(value, bool):
            return 0.5 if value else -0.5
        if isinstance(value, (int, float)):
            if math.isclose(value, 0.0, rel_tol=1e-3, abs_tol=1e-3):
                return -0.5
            if math.isclose(value, 1.0, rel_tol=1e-3, abs_tol=1e-3):
                return 0.5
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in {"y", "yes", "pass", "true", "+"}:
                return 0.5
            if lowered in {"n", "no", "fail", "false", "-"}:
                return -0.5
    return float(value)


def _normalized_feature_value(feature: str, value: float | None) -> float | None:
    """将原始特征值映射到 0-1 区间；缺失值返回 None。"""

    prepared = _prepare_feature_value(feature, value)
    if prepared is None:
        return None
    window = FEATURE_WINDOWS.get(feature)
    if window is None:
        return None
    return _window_score(prepared, window)


def compute_desirability_score(ingredient: Ingredient) -> float:
    """计算单体可用性分（几何平均）。"""

    scores: list[float] = []
    for feature, attr in ADME_FEATURE_ATTRS.items():
        value = getattr(ingredient, attr)
        normalized = _normalized_feature_value(feature, value)
        if normalized is not None:
            scores.append(max(normalized, _SAFE_EPS))
    if not scores:
        return 0.0
    product = 1.0
    for score in scores:
        product *= score
    # 采用几何平均放大小值影响，使任一指标严重偏离时整体得分下降。
    geo_mean = product ** (1 / len(scores))
    coverage = len(scores) / max(_TOTAL_FEATURES, 1)
    coverage_floor = max(0.0, min(1.0, _get_aci_settings().get("coverage_floor", 0.0)))
    return geo_mean * max(coverage, coverage_floor)


def compute_complementarity_score(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient]
) -> float:
    """根据特征多样性计算组合互补性。"""

    feature_diversities: list[float] = []
    variance_cap = max(
        _SAFE_EPS,
        float(_get_aci_settings().get("complement_variance_max", 0.25)),
    )
    for feature, attr in ADME_FEATURE_ATTRS.items():
        values: list[float] = []
        weights: list[float] = []
        for select, proportion, ingredient in zip(
            candidate.iter_selects(), candidate.proportions, ingredients
        ):
            if not select:
                continue
            normalized = _normalized_feature_value(feature, getattr(ingredient, attr))
            if normalized is None:
                continue
            values.append(normalized)
            weights.append(proportion)
        if len(values) <= 1:
            continue
        import numpy as np  # 惰性导入，避免默认运行时依赖 numpy
        weights_arr = np.array(weights, dtype=float)
        total = weights_arr.sum()
        if total <= 0:
            continue
        weights_arr /= total
        values_arr = np.array(values, dtype=float)
        mean = float(np.dot(weights_arr, values_arr))
        variance = float(np.dot(weights_arr, (values_arr - mean) ** 2))
        feature_diversities.append(min(1.0, variance / variance_cap))
    if not feature_diversities:
        return 0.0
    return float(sum(feature_diversities) / len(feature_diversities))


def compute_aci_score(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    lambda_weight: float | None = None,
) -> float:
    """计算基于 ADME 特征的 ACI 指标。"""

    desirabilities: list[tuple[float, float]] = []
    for select, proportion, ingredient in zip(
        candidate.iter_selects(), candidate.proportions, ingredients
    ):
        if not select:
            continue
        desirability = compute_desirability_score(ingredient)
        desirabilities.append((proportion, desirability))
    if not desirabilities:
        return 0.0
    total_prop = sum(weight for weight, _ in desirabilities)
    if total_prop <= 0:
        return 0.0
    weighted_mean = (
        sum(weight * desirability for weight, desirability in desirabilities)
        / total_prop
    )
    complementarity = compute_complementarity_score(candidate, ingredients)
    if lambda_weight is None:
        lambda_weight = _get_aci_settings().get("lambda_weight", ACI_LAMBDA)
    lambda_weight = max(0.0, min(1.0, lambda_weight))
    return (1 - lambda_weight) * weighted_mean + lambda_weight * complementarity


def compute_hepatotoxicity_score(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient]
) -> float:
    """计算组合的加权肝毒性得分。

    Args:
        candidate: 已归一化组合。
        ingredients: 成分数据集。

    Returns:
        配比加权后的 hepatotoxicity_score。
    """
    score = 0.0
    for select, proportion, ingredient in zip(
        candidate.iter_selects(), candidate.proportions, ingredients
    ):
        if select:
            score += ingredient.hepatotoxicity_score * proportion
    return score


def diversity_penalty(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient]
) -> float:
    """若单一成分占比 ≥80%，返回 0.1 惩罚，否则 0。

    Args:
        candidate: 当前组合。
        ingredients: 成分信息（占位，保留与 evaluate_metrics 的签名一致）。

    Returns:
        惩罚值（0.1 或 0）。
    """
    max_ratio = 0.0
    for select, proportion in zip(candidate.iter_selects(), candidate.proportions):
        if select:
            max_ratio = max(max_ratio, proportion)
    return 0.1 if max_ratio >= 0.8 else 0.0


def evaluate_metrics(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient]
) -> Tuple[float, float, float, float]:
    """返回 (ACI(含惩罚), 肝毒性, 惩罚, 纯 ACI)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (ACI, 肝毒性, 惩罚, ACI_raw)
    """
    # 先归一化配比，再统一计算各项指标，确保多算法输出可比较。
    normalized_candidate = candidate.with_normalized()
    aci_raw = compute_aci_score(normalized_candidate, ingredients)
    toxicity = compute_hepatotoxicity_score(normalized_candidate, ingredients)
    penalty = diversity_penalty(normalized_candidate, ingredients)
    aci_with_penalty = aci_raw * (1 - penalty)
    return aci_with_penalty, toxicity, penalty, aci_raw


def single_objective_score(
    candidate: CandidateSolution,
    ingredients: Sequence[Ingredient],
    toxicity_weight: float = 1.0,
) -> float:
    """单目标评价：ACI - toxicity_weight * 肝毒性。

    Args:
        candidate: 当前候选。
        ingredients: 成分列表。
        toxicity_weight: 毒性在适应度中的相对权重。

    Returns:
        复方综合得分（越高越优）。
    """
    aci, toxicity, _, _ = evaluate_metrics(candidate, ingredients)
    return aci - toxicity_weight * toxicity


def multi_objective_score(
    candidate: CandidateSolution, ingredients: Sequence[Ingredient]
) -> Tuple[float, float]:
    """返回 (ACI, 肝毒性)。

    Args:
        candidate: 当前组合。
        ingredients: 成分列表。

    Returns:
        (ACI, 肝毒性)。
    """
    aci, toxicity, _, _ = evaluate_metrics(candidate, ingredients)
    return aci, toxicity


class CombinationProblem(_ElementwiseProblem):
    """pymoo 兼容问题定义，用于描述复方的选择+配比变量与 ACI/毒性目标。

    Args:
        ingredients: 候选成分列表。
        mode: \"single\" 表示单目标（ACI - 毒性），\"multi\" 表示多目标。
        toxicity_weight: 仅在单目标模式下生效，表示肝毒性惩罚权重。
    """

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        mode: str = "single",
        toxicity_weight: float = 1.0,
    ):
        n_var = len(ingredients) * 2
        n_obj = 1 if mode == "single" else 2
        # 变量下界/上界均为 [0,1]，以便直接套用向量编码。
        super().__init__(n_var=n_var, n_obj=n_obj, xl=0.0, xu=1.0)
        self.ingredients = ingredients
        self.mode = mode
        self.toxicity_weight = toxicity_weight

    def _evaluate(self, x, out, *args, **kwargs):
        """在 pymoo 运行过程中评估单个向量。

        Args:
            x: 当前向量（shape: n_var）。
            out: 输出 dict，pymoo 要求填入 \"F\"。
        """
        import numpy as np  # 惰性导入
        candidate = vector_to_candidate(x, self.ingredients)
        aci, toxicity = multi_objective_score(candidate, self.ingredients)
        if self.mode == "single":
            # pymoo 默认最小化，故取负值以最大化 (ACI - 权重*毒性)。
            out["F"] = np.array([-(aci - self.toxicity_weight * toxicity)])
        else:
            out["F"] = np.array([-aci, toxicity])


class PymooSingleObjectiveGA:
    """使用 pymoo GA 处理单目标优化（ACI - 毒性）。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        toxicity_weight: float | None = None,
        generations: int | None = None,
        population_size: int | None = None,
    ):
        """配置 GA 所需组件；默认参数来自 config/algorithms.json。"""
        from pymoo.algorithms.soo.nonconvex.ga import GA as PymooGAAlg  # type: ignore
        cfg = _load_algorithm_config()["ga"]
        self.toxicity_weight = (
            toxicity_weight if toxicity_weight is not None else cfg["toxicity_weight"]
        )
        self.generations = (
            generations if generations is not None else cfg["generations"]
        )
        pop = population_size if population_size is not None else cfg["population_size"]
        self.problem = CombinationProblem(
            ingredients, mode="single", toxicity_weight=self.toxicity_weight
        )
        self.algorithm = PymooGAAlg(pop_size=pop, eliminate_duplicates=True)

    def run(self) -> CandidateSolution:
        """执行单目标 GA 并返回最优 candidate。"""
        from pymoo.optimize import minimize  # type: ignore
        result = minimize(
            self.problem,
            self.algorithm,
            termination=("n_gen", self.generations),
            verbose=False,
        )
        return vector_to_candidate(result.X, self.problem.ingredients)


class PymooNSGAII:
    """基于 pymoo NSGA-II 的多目标优化器。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        generations: int | None = None,
        population_size: int | None = None,
    ):
        """初始化 NSGA-II；默认参数读取配置文件。"""
        from pymoo.algorithms.moo.nsga2 import NSGA2  # type: ignore
        cfg = _load_algorithm_config()["nsga2"]
        self.generations = (
            generations if generations is not None else cfg["generations"]
        )
        pop = population_size if population_size is not None else cfg["population_size"]
        self.problem = CombinationProblem(ingredients, mode="multi")
        self.algorithm = NSGA2(pop_size=pop)

    def run(self) -> Tuple[List[CandidateSolution], "np.ndarray"]:
        """执行 NSGA-II 并返回 candidate 与目标值矩阵。

        Returns:
            solutions: 输入向量转化的 CandidateSolution 列表。
            np.ndarray: 每个 solution 的目标值（ACI 取负 + 肝毒性）。
        """
        from pymoo.optimize import minimize  # type: ignore
        result = minimize(
            self.problem,
            self.algorithm,
            termination=("n_gen", self.generations),
            verbose=False,
        )
        solutions = [vector_to_candidate(x, self.problem.ingredients) for x in result.X]
        return solutions, result.F


class PySwarmsPSO:
    """借助 pyswarms 实现的单目标粒子群优化。"""

    def __init__(
        self,
        ingredients: Sequence[Ingredient],
        iterations: int | None = None,
        swarm_size: int | None = None,
        options: dict | None = None,
    ):
        """配置 PSO 超参数，默认来自配置文件。"""
        cfg = _load_algorithm_config()["pso"]
        self.ingredients = ingredients
        self.iterations = iterations if iterations is not None else cfg["iterations"]
        self.swarm_size = swarm_size if swarm_size is not None else cfg["swarm_size"]
        self.options = options if options is not None else cfg["options"]
        self.dimensions = len(ingredients) * 2

    def _fitness(self, swarm: "np.ndarray") -> "np.ndarray":
        """pyswarms 回调：计算每个粒子的目标值（越小越优）。"""
        import numpy as np  # 惰性导入
        scores = []
        for particle in swarm:
            candidate = vector_to_candidate(particle, self.ingredients)
            score = single_objective_score(candidate, self.ingredients)
            scores.append(-score)
        return np.array(scores)

    def run(self) -> CandidateSolution:
        """执行 PSO 并返回最优 candidate。"""
        import numpy as np  # 惰性导入
        from pyswarms.single.global_best import GlobalBestPSO  # type: ignore
        optimizer = GlobalBestPSO(
            n_particles=self.swarm_size,
            dimensions=self.dimensions,
            options=self.options,
            bounds=(np.zeros(self.dimensions), np.ones(self.dimensions)),
        )
        _, best_pos = optimizer.optimize(
            self._fitness, iters=self.iterations, verbose=False
        )
        return vector_to_candidate(best_pos, self.ingredients)


def select_candidate_by_weighted_sum(
    pareto_candidates: Sequence[CandidateSolution],
    ingredients: Sequence[Ingredient],
    toxicity_weight: float,
) -> Tuple[CandidateSolution, dict]:
    """从 Pareto 集合中按权重（ACI-毒性）挑选最优解。

    Returns:
        (最佳 candidate, 指标字典)。
    """

    if not pareto_candidates:
        raise ValueError("Pareto 集合为空，无法做权重挑选。")
    best = None
    best_score = -math.inf
    for candidate in pareto_candidates:
        score = single_objective_score(candidate, ingredients, toxicity_weight)
        if score > best_score:
            best = candidate
            best_score = score
    assert best is not None
    aci_pen, tox, penalty, aci_raw = evaluate_metrics(best, ingredients)
    return best, {
        "score": best_score,
        "aci_penalized": aci_pen,
        "aci_raw": aci_raw,
        "toxicity": tox,
        "penalty": penalty,
    }
