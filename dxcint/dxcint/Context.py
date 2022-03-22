import dxpy
import inspect


class ContextError(Exception):
    """
    Class to handle Context errors
    """


class Context(object):
    def __init__(
            self,
            project: str
    ):
        self._project_id = self._resolve_project(project=project)

    def __setattr__(self, *args):
        if inspect.stack()[1][3] == "__init__":
            object.__setattr__(self, *args)
        else:
            raise ContextError("Context class is immutable")

    @property
    def project_id(self):
        return self._project_id

    @staticmethod
    def _resolve_project(project: str) -> str:
        # First, see if the project is a project-id.
        try:
            project = dxpy.DXProject(project)
            return project.get_id()
        except dxpy.DXError:
            found_projects = list(
                dxpy.find_projects(name=project, level="VIEW")
            )
            if len(found_projects) == 0:
                raise ContextError(
                    f"Context._resolve_project(): Could not find project `{project}`."
                )
            elif len(found_projects) > 1:
                raise ContextError(
                    f"Context._resolve_project(): found multiple projects with name `{project}`"
                )
            else:
                return found_projects[0].get("id")
