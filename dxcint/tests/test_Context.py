import dxpy
from dxcint.Context import Context


def test__resolve_project():
    context = Context(project="dxCompiler_playground")
    assert context.project_id == list(dxpy.find_projects(name="dxCompiler_playground", level="VIEW"))[0].get("id")
