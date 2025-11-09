import ast
from pathlib import Path


EXCLUDE = {".venv", ".pytest_cache", ".serena", "__pycache__"}


def _should_skip(path: Path) -> bool:
    return any(part in EXCLUDE for part in path.parts)


def _literal_info(node: ast.AST):
    if isinstance(node, ast.Constant):
        if node.value is None:
            return None
        return repr(node.value)
    if isinstance(node, (ast.List, ast.Tuple, ast.Set, ast.Dict)):
        try:
            # Python 3.11+ 可用
            return ast.unparse(node)  # type: ignore[attr-defined]
        except Exception:
            return node.__class__.__name__
    return None


def test_no_is_with_literal():
    root = Path(__file__).resolve().parent.parent
    issues = []
    for path in root.rglob("*.py"):
        rel = path.relative_to(root)
        if _should_skip(rel):
            continue
        text = path.read_text(encoding="utf-8")
        tree = ast.parse(text, filename=str(rel))
        lines = text.splitlines()

        class V(ast.NodeVisitor):
            def visit_Compare(self, node: ast.Compare):
                operands = [node.left, *node.comparators]
                for i, op in enumerate(node.ops):
                    if isinstance(op, (ast.Is, ast.IsNot)):
                        l = operands[i]
                        r = operands[i + 1]
                        li = _literal_info(l)
                        ri = _literal_info(r)
                        if li or ri:
                            issues.append(
                                (
                                    str(rel),
                                    node.lineno,
                                    "is" if isinstance(op, ast.Is) else "is not",
                                    li,
                                    ri,
                                    lines[node.lineno - 1].strip(),
                                )
                            )
                self.generic_visit(node)

        V().visit(tree)

    assert not issues, "发现使用 'is/is not' 与字面量比较的问题:\n" + "\n".join(
        f"{rel}:{ln}:{op} left={li!r} right={ri!r} | {line}"
        for rel, ln, op, li, ri, line in issues
    )

