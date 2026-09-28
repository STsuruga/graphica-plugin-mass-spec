# tests/conftest.py
"""テスト共通のフィクスチャ。

Graphica 本体(pip install "graphica-plot>=2.0,<3")が無い環境では、本体に依存しない
計算のテストだけが走り、配線・読み込みのテストは skip される。
"""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

try:
    from PySide6.QtWidgets import QApplication
except ImportError:  # pragma: no cover - PySide6が無い環境
    QApplication = None

try:
    import graphica.plugin.testing  # noqa: F401
    GRAPHICA_AVAILABLE = True
except ImportError:
    GRAPHICA_AVAILABLE = False


requires_graphica = pytest.mark.skipif(
    not GRAPHICA_AVAILABLE,
    reason="Graphica 本体が未インストールです(pip install \"graphica-plot>=2.0,<3\")",
)


@pytest.fixture(scope="session", autouse=True)
def qapp():
    """Graphica 本体と同じく、セッション全体で QApplication を1つだけ用意する。"""
    if QApplication is None:
        yield None
        return
    app = QApplication.instance() or QApplication([])
    yield app
