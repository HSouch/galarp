import nox


@nox.session
def tests(session):
    """ Checking style using ruff"""
    session.install("-r", "requirements.txt")

    session.run("pytest")

@nox.session
def lint(session):
    """ Checking style using flake8"""
    session.install("ruff")
    session.run("ruff", "check", "galarp")