import nox


@nox.session
def lint(session):
    """ Checking style using ruff"""
    
    session.install()
    session.run("pytest")

    #session.install("ruff")
    #session.run("ruff check", "galarp")
