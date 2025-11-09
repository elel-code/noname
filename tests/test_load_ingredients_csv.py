import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from algorithms import Ingredient, load_ingredients_csv


def test_load_ingredients_csv_minimal(tmp_path: Path):
    csv_text = (
        "name,smiles,ob,dl,mw,alogp,h_don,h_acc,caco2,bbb,fasa,hl\n"
        "黄芩素,OC1=CC=C2C(=C1)C(=O)C=C(O)C(=O)O2,35,0.18,270.2,2.1,3,5,0.82,-0.3,0.33,6.5\n"
        "川芎嗪,CN1C=NC(=CN=C1N(C)C)N(C)C,32,0.09,136.2,0.3,1,2,1.1,1,0.28,3.5\n"
    )
    p = tmp_path / "ingredients.csv"
    p.write_text(csv_text, encoding="utf-8")

    items = load_ingredients_csv(p)

    assert isinstance(items, list) and len(items) == 2
    assert all(isinstance(x, Ingredient) for x in items)
    assert items[0].name == "黄芩素"
    # 未提供 hepatotoxicity_score，应由 __post_init__ 自动推断填充
    assert items[0].hepatotoxicity_score is not None
