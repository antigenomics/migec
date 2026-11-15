"""Interactive branding must leave redirected output and machine-readable stdout clean."""
from migec.cli import main


def test_redirected_stderr_is_silent(monkeypatch, capsys):
    monkeypatch.setattr("sys.stderr.isatty", lambda: False)
    main()
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""


def test_terminal_splash_uses_stderr(monkeypatch, capsys):
    monkeypatch.setattr("sys.stderr.isatty", lambda: True)
    monkeypatch.setenv("TERM", "xterm-256color")
    main()
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "migec" in captured.err
    assert "checkout" in captured.err
