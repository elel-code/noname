from __future__ import annotations

"""应用编排模块：组织单/多目标流程与输出。

- 惰性导入算法与可视化，降低默认运行成本；
- 提供统一的文本摘要、Top-K、导出与图表保存；
- 面向 CLI（main.py）进行功能编排与复用。
"""

import json
import sys
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:  # 仅类型检查与补全，运行时不引入依赖
    from algorithms import CandidateSolution, Ingredient
    import numpy as np


def sample_ingredients() -> list[Ingredient]:
    """返回示例成分列表，用于算法演示。"""

    # 延迟导入，避免默认模式下的依赖报错
    from algorithms import Ingredient

    return [
        Ingredient(
            "黄芩素",
            smiles="OC1=CC=C2C(=C1)C(=O)C=C(O)C(=O)O2",
            hepatotoxicity_score=0.3,
            ob=35.0,
            dl=0.18,
            mw=270.2,
            alogp=2.1,
            h_don=3,
            h_acc=5,
            caco2=0.82,
            bbb=-0.3,
            fasa=0.33,
            hl=6.5,
        ),
        Ingredient(
            "川芎嗪",
            smiles="CN1C=NC(=CN=C1N(C)C)N(C)C",
            hepatotoxicity_score=0.2,
            ob=32.0,
            dl=0.09,
            mw=136.2,
            alogp=0.3,
            h_don=1,
            h_acc=2,
            caco2=1.1,
            bbb=0.2,
            fasa=0.28,
            hl=3.5,
        ),
        Ingredient(
            "黄连碱",
            smiles="COC1=CC=CC=C1O",
            ob=26.0,
            dl=0.12,
            mw=336.4,
            alogp=1.5,
            h_don=0,
            h_acc=4,
            caco2=0.4,
            bbb=-0.7,
            fasa=0.48,
            hl=8.0,
        ),
        Ingredient(
            "丹参酮",
            smiles="CC1=C2C(=O)C=C(C2=CC3=CC=CC=C31)C",
            hepatotoxicity_score=0.15,
            ob=45.0,
            dl=0.26,
            mw=294.3,
            alogp=3.1,
            h_don=0,
            h_acc=4,
            caco2=1.2,
            bbb=0.5,
            fasa=0.22,
            hl=10.5,
        ),
        Ingredient(
            "甘草酸",
            smiles="CC1=CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C(=O)O",
            hepatotoxicity_score=0.4,
            ob=28.0,
            dl=0.11,
            mw=822.0,
            alogp=1.2,
            h_don=5,
            h_acc=16,
            caco2=0.2,
            bbb=-1.2,
            fasa=0.62,
            hl=4.0,
        ),
        Ingredient(
            "白芍苷",
            smiles="COC1=CC=C2C(=O)C3=C(C=C(C(=C3O2)O)O)O1",
            hepatotoxicity_score=0.18,
            ob=38.0,
            dl=0.14,
            mw=480.4,
            alogp=-1.1,
            h_don=6,
            h_acc=11,
            caco2=0.35,
            bbb=-0.9,
            fasa=0.58,
            hl=5.5,
        ),
    ]


def describe_candidate(label: str, candidate: CandidateSolution, ingredients: list[Ingredient]) -> None:
    """打印候选结果的指标与具体成分选取。"""

    # 延迟导入以避免默认模式触发依赖
    from algorithms import evaluate_metrics

    aci_penalized, toxicity, penalty, aci_raw = evaluate_metrics(candidate, ingredients)
    normalized = candidate.with_normalized()
    mix = []
    for idx in normalized.iter_selected_indices():
        percentage = normalized.proportions[idx] * 100
        mix.append(f"{ingredients[idx].name}({percentage:.1f}%)")
    mix_repr = "[" + ", ".join(mix) + "]" if mix else "[]"
    # 对外展示时使用 penalized ACI，避免占比过度集中。
    print(
        f"{label}：ACI_pen {aci_penalized:.3f} (raw {aci_raw:.3f})，"
        f"肝毒性 {toxicity:.3f}，惩罚 {penalty:.2f}，配伍 {mix_repr}"
    )


def _format_mix(candidate: CandidateSolution, ingredients: list[Ingredient], *, percent_digits: int = 1) -> str:
    """返回候选解的配伍配比字符串，如 [黄芩素(30.0%), 丹参酮(70.0%)]。"""
    normalized = candidate.with_normalized()
    parts: list[str] = []
    for idx in normalized.iter_selected_indices():
        pct = normalized.proportions[idx] * 100
        parts.append(f"{ingredients[idx].name}({pct:.{percent_digits}f}%)")
    return "[" + ", ".join(parts) + "]" if parts else "[]"


def _extract_mix(candidate: CandidateSolution, ingredients: list[Ingredient]) -> list[dict[str, float]]:
    """结构化提取配伍配比，返回 [{name, percent}, ...]，percent 为 0..100。"""
    normalized = candidate.with_normalized()
    items: list[dict[str, float]] = []
    for idx in normalized.iter_selected_indices():
        pct = float(normalized.proportions[idx] * 100.0)
        items.append({"name": ingredients[idx].name, "percent": pct})
    return items


def run_single_objective_algorithms(
    ingredients: list[Ingredient],
) -> tuple[list[CandidateSolution], float, dict]:
    """执行 pymoo GA 与 pyswarms PSO，并返回两个最优解与收敛历史。"""

    # 延迟导入
    from algorithms import PymooSingleObjectiveGA, PySwarmsPSO

    ga = PymooSingleObjectiveGA(ingredients)
    pso = PySwarmsPSO(ingredients)
    best_ga = ga.run()
    best_pso = pso.run()
    info = {
        "ga_hist": getattr(ga, "convergence_", None),
        "pso_hist": getattr(pso, "convergence_", None),
    }
    return [best_ga, best_pso], ga.toxicity_weight, info


def _subsample_solutions_by_toxicity(
    solutions: list[CandidateSolution], metrics: "np.ndarray", k: int
):
    """按肝毒性排序等间距抽取最多 k 个候选，控制可视化与打印体量。

    说明：
    - metrics[:,1] 为毒性（越低越好），我们在全排序后取线性等间距下标，保留两端；
    - 仅在候选数量明显多于输出上限时触发。
    """
    import numpy as np

    if k is None or k <= 0:
        return solutions, metrics
    n = len(solutions)
    if k >= n:
        return solutions, metrics
    # 毒性升序排序后按等间距抽样
    order = np.argsort(metrics[:, 1])
    idxs = np.linspace(0, n - 1, num=k).round().astype(int)
    chosen = order[idxs]
    chosen = np.unique(chosen)
    return [solutions[i] for i in chosen], metrics[chosen]


def run_nsga(
    ingredients: list[Ingredient],
    toxicity_weight: float,
    highlight_points: list[tuple[str, float, float]] | None = None,
    *,
    max_candidates: int | None = None,
    precision: int = 3,
    save_metrics: str | None = None,
    viz_enabled: bool = True,
    viz_output_dir: str | None = None,
    viz_dpi: int | None = None,
    viz_colormap: str | None = None,
    viz_figs: list[str] | None = None,
    viz_formats: list[str] | None = None,
    viz_annotate: bool | None = None,
    viz_responsive: bool | None = None,
    viz_title_prefix: str | None = None,
    summary_top_k: int | None = 3,
    export_json: str | None = None,
) -> None:
    """运行 NSGA-II 并打印非支配集候选，生成可视化与汇总。

    参数：
    - ingredients: 成分列表（算法搜索空间由其长度决定）；
    - toxicity_weight: 最终从 Pareto 集中按 ACI-权重*毒性进行加权挑选的权重；
    - highlight_points: 可在帕累托图中高亮的锚点（名称、毒性、ACI_raw）；
    - max_candidates: 候选上限（过多时抽样，便于阅读与作图）；
    - precision: 文本输出保留的小数位；
    - save_metrics: 若提供则把候选指标写入 CSV；
    - viz_*: 控制是否出图、目录、DPI、配色、图型集合、保存格式、注释、自适应、标题前缀；
    - summary_top_k: 打印 Top-K（按 ACI 与最低毒性）；
    - export_json: 导出机器可读汇总 JSON 的路径。
    """

    # 延迟导入
    from algorithms import PymooNSGAII, evaluate_metrics, select_candidate_by_weighted_sum
    from visualization import render_nsga_visualizations

    nsga = PymooNSGAII(ingredients)
    solutions, metrics = nsga.run()
    if max_candidates:
        solutions, metrics = _subsample_solutions_by_toxicity(solutions, metrics, max_candidates)

    print("\nNSGA-II 非支配集候选：")
    for idx, (candidate, metric) in enumerate(zip(solutions, metrics), 1):
        aci = -metric[0]
        toxicity = metric[1]
        penalty = evaluate_metrics(candidate, ingredients)[2]
        mix_repr = _format_mix(candidate, ingredients, percent_digits=1)
        print(
            f"  {idx}. ACI {aci:.{precision}f}，肝毒性 {toxicity:.{precision}f}，惩罚 {penalty:.2f}，配伍 {mix_repr}"
        )

    if save_metrics:
        try:
            import csv

            path = Path(save_metrics)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", newline="", encoding="utf-8") as fp:
                writer = csv.writer(fp)
                # 增加 mix 列，便于批处理查看配伍配比
                writer.writerow(["# ACI(+)", "toxicity", "penalty", "mix"])
                for cand, m in zip(solutions, metrics):
                    aci = -m[0]
                    tox = m[1]
                    pen = evaluate_metrics(cand, ingredients)[2]
                    mix_repr = _format_mix(cand, ingredients, percent_digits=1)
                    writer.writerow([f"{aci:.{precision}f}", f"{tox:.{precision}f}", f"{pen:.2f}", mix_repr])
            print(f"\n候选指标已保存：{path}")
        except Exception as e:  # noqa: BLE001
            print(f"保存指标失败：{e}")

    best_candidate, best_metrics = select_candidate_by_weighted_sum(solutions, ingredients, toxicity_weight)
    print(
        f"\n按 toxicity_weight={toxicity_weight:.2f} 的加权结果："
        f"ACI_pen {best_metrics['aci_penalized']:.{precision}f} (raw {best_metrics['aci_raw']:.{precision}f})，"
        f"肝毒性 {best_metrics['toxicity']:.{precision}f}，惩罚 {best_metrics['penalty']:.2f}"
    )
    describe_candidate("NSGA-权重选择", best_candidate, ingredients)

    # 汇总统计与 Top-K
    try:
        import numpy as _np

        aci_arr = _np.array([-m[0] for m in metrics], dtype=float)
        tox_arr = _np.array([m[1] for m in metrics], dtype=float)
        print("\n统计汇总：")
        print(
            f"  候选数 {len(solutions)} | ACI[min {aci_arr.min():.{precision}f}, max {aci_arr.max():.{precision}f}, mean {aci_arr.mean():.{precision}f}]"
        )
        print(
            f"  肝毒性[min {tox_arr.min():.{precision}f}, max {tox_arr.max():.{precision}f}, mean {tox_arr.mean():.{precision}f}]"
        )
        k = summary_top_k or 0
        if k > 0:
            top_aci_idx = _np.argsort(-aci_arr)[:k]
            top_safe_idx = _np.argsort(tox_arr)[:k]
            print("\nTop-K（按 ACI）：")
            for rank, i in enumerate(top_aci_idx, 1):
                mix_repr = _format_mix(solutions[i], ingredients, percent_digits=1)
                print(
                    f"  {rank}. idx {i+1} → ACI {aci_arr[i]:.{precision}f}，毒性 {tox_arr[i]:.{precision}f}，配伍 {mix_repr}"
                )
            print("Top-K（按最低毒性）：")
            for rank, i in enumerate(top_safe_idx, 1):
                mix_repr = _format_mix(solutions[i], ingredients, percent_digits=1)
                print(
                    f"  {rank}. idx {i+1} → 毒性 {tox_arr[i]:.{precision}f}，ACI {aci_arr[i]:.{precision}f}，配伍 {mix_repr}"
                )
    except Exception:
        pass

    # 高亮点（用于图中标签）
    highlight_points = list(highlight_points or [])
    highlight_points.append(("NSGA-权重", best_metrics["toxicity"], best_metrics["aci_raw"]))

    if viz_enabled:
        artifacts = render_nsga_visualizations(
            best_candidate,
            ingredients,
            metrics,
            highlight_points,
            output_dir=Path(viz_output_dir) if viz_output_dir else None,
            dpi=viz_dpi or 200,
            colormap=viz_colormap or "viridis",
            figs=tuple(viz_figs) if viz_figs else None,
            formats=tuple(viz_formats) if viz_formats else None,
            annotate=True if viz_annotate is None else bool(viz_annotate),
            responsive=True if viz_responsive is None else bool(viz_responsive),
            title_prefix=viz_title_prefix,
        )
        if artifacts:
            joined = "、".join(str(p) for p in artifacts.values())
            print(f"\n可视化已保存：{joined}")

    # 导出 JSON 汇总（含逐候选配伍配比与最佳权重解的配伍）
    if export_json:
        try:
            aci_list = [-float(m[0]) for m in metrics]
            tox_list = [float(m[1]) for m in metrics]
            sols: list[dict[str, Any]] = []
            for cand, m in zip(solutions, metrics):
                aci = -float(m[0])
                tox = float(m[1])
                pen = float(evaluate_metrics(cand, ingredients)[2])
                sols.append(
                    {
                        "aci": aci,
                        "toxicity": tox,
                        "penalty": pen,
                        "mix": _extract_mix(cand, ingredients),
                    }
                )
            payload: dict[str, Any] = {
                "count": len(aci_list),
                "aci": aci_list,
                "toxicity": tox_list,
                "solutions": sols,
                "best_weighted": {
                    **best_metrics,
                    "mix": _extract_mix(best_candidate, ingredients),
                },
            }
            out_path = Path(export_json)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
            print(f"汇总 JSON 已导出：{out_path}")
        except Exception as e:  # noqa: BLE001
            print(f"导出 JSON 失败：{e}")


def run_demo(ingredients: list["Ingredient"] | None = None, *, args=None) -> int:
    """运行演示（包含算法与可视化）。"""

    try:
        from visualization import configure_matplotlib_fonts
    except Exception:
        print(
            "缺少可视化或相关依赖，无法运行 demo。\n"
            "请在可写环境下安装依赖后重试（例如：pymoo、pyswarms、matplotlib 等）",
            file=sys.stderr,
        )
        return 2

    if not (args and args.no_fig):
        preferred = ([args.font] if (args and args.font) else None)
        fallback = (Path(args.font_dir) if (args and args.font_dir) else None)
        configure_matplotlib_fonts(preferred_fonts=preferred, fallback_dir=fallback)

    ingredients = ingredients or sample_ingredients()
    single_objective_results, toxicity_weight, so_info = run_single_objective_algorithms(ingredients)
    print("单目标算法结果：")
    for label, candidate in zip(("遗传算法", "粒子群算法"), single_objective_results):
        describe_candidate(label, candidate, ingredients)

    # 延迟导入评估函数，避免默认模式触发依赖
    from algorithms import evaluate_metrics

    highlight_points: list[tuple[str, float, float]] = []
    for label, candidate in zip(("GA 最优", "PSO 最优"), single_objective_results):
        _, toxicity, _, aci_raw = evaluate_metrics(candidate, ingredients)
        highlight_points.append((label, toxicity, aci_raw))

    # 读取可视化默认配置
    try:
        from algorithms import _load_algorithm_config

        cfg = _load_algorithm_config()
        viz_cfg = cfg.get("viz", {}) if isinstance(cfg, dict) else {}
    except Exception:
        cfg = {}
        viz_cfg = {}

    # 单目标可视化（各自独立图）
    if not (args and args.no_fig):
        try:
            from visualization import render_single_objective_visualizations

            so_artifacts = render_single_objective_visualizations(
                best_ga=single_objective_results[0],
                best_pso=single_objective_results[1],
                ingredients=ingredients,
                ga_convergence=so_info.get("ga_hist"),
                pso_convergence=so_info.get("pso_hist"),
                output_dir=Path(args.output_dir) if (args and args.output_dir) else None,
                dpi=(args.dpi if args and args.dpi else 200),
                formats=(
                    [s.strip().lower() for s in args.img_formats.split(",")]
                    if (args and getattr(args, "img_formats", None))
                    else None
                ),
                title_prefix=(args.title_prefix if args and args.title_prefix else None),
            )
            if so_artifacts:
                joined = "、".join(str(p) for p in so_artifacts.values())
                print(f"单目标可视化已保存：{joined}")
        except Exception:
            pass

    run_nsga(
        ingredients,
        toxicity_weight,
        highlight_points,
        max_candidates=(
            args.max_candidates
            if args and args.max_candidates is not None
            else viz_cfg.get("max_candidates") or cfg.get("nsga2", {}).get("max_candidates")
            if "cfg" in locals()
            else None
        ),
        precision=(args.precision if args else 3),
        save_metrics=(args.save_metrics if args else None),
        viz_enabled=(not (args and args.no_fig)) and bool(viz_cfg.get("enabled", True)),
        viz_output_dir=(args.output_dir if args and args.output_dir else viz_cfg.get("output_dir")),
        viz_dpi=(args.dpi if args and args.dpi else viz_cfg.get("dpi")),
        viz_colormap=(args.colormap if args and args.colormap else viz_cfg.get("colormap")),
        viz_figs=(
            [s.strip().lower() for s in args.figs.split(",")]
            if (args and getattr(args, "figs", None))
            else (list(viz_cfg.get("figs")) if isinstance(viz_cfg.get("figs"), (list, tuple)) else None)
        ),
        viz_formats=(
            [s.strip().lower() for s in args.img_formats.split(",")]
            if (args and getattr(args, "img_formats", None))
            else (
                list(viz_cfg.get("img_formats")) if isinstance(viz_cfg.get("img_formats"), (list, tuple)) else None
            )
        ),
        viz_annotate=(False if (args and getattr(args, "no_annotate", False)) else viz_cfg.get("annotate", True)),
        viz_responsive=(
            True
            if (args and getattr(args, "responsive", False))
            else False
            if (args and getattr(args, "no_responsive", False))
            else viz_cfg.get("responsive", True)
        ),
        viz_title_prefix=(args.title_prefix if args and args.title_prefix else viz_cfg.get("title_prefix")),
        summary_top_k=(args.summary_top_k if args and getattr(args, "summary_top_k", None) is not None else 3),
        export_json=(args.export_json if args and getattr(args, "export_json", None) else None),
    )
    return 0
