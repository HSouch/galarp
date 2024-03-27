import nox


@nox.session
def tests(session):
    """ Checking style using ruff"""
    session.install("-r", "requirements.txt")

    session.run("pytest")

    #session.install("ruff")
    #session.run("ruff check", "galarp")
