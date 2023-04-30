import dxpy
import itertools
from typing import Optional, Dict, Generator, Set, Iterable, Union

from dxcint.RegisteredTest import RegisteredTest


class JobCollectorMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _collect(self) -> Set[str]:
        """
        Collects unique job IDs spawned by the provided execution.
        Returns: Set(str). Job IDs spawned by the execution provided upon class instantiation.
        """
        child_executions = dxpy.find_executions(
            project=self.context.project_id, **self._prep_find_kwargs(self.job_id)
        )
        return set(self._recursive_collect(child_executions))

    def _recursive_collect(
        self, execution_ids: Union[Generator, Iterable]
    ) -> Generator:
        for this_execution in execution_ids:
            this_execution_id = this_execution.get("id", "")
            if this_execution_id.startswith("analysis-"):
                unpacked_new_ids = dxpy.find_executions(
                    **self._prep_find_kwargs(this_execution_id)
                )
                yield from self._recursive_collect(unpacked_new_ids)
            elif this_execution_id.startswith("job-"):
                yield this_execution_id
                unpacked_new_ids = dxpy.find_executions(
                    **self._prep_find_kwargs(this_execution_id)
                )
                peeked_generator = self._peek_generator(unpacked_new_ids)
                if peeked_generator:
                    yield from self._recursive_collect(peeked_generator)
            else:
                raise ValueError(
                    f"Unknown ID type should be job or analysis. Provided {this_execution_id}"
                )

    @staticmethod
    def _peek_generator(gen: Generator) -> Optional[Iterable]:
        try:
            first = next(gen)
        except StopIteration:
            return None
        return itertools.chain([first], gen)

    @staticmethod
    def _prep_find_kwargs(execution: str) -> Dict:
        if execution.startswith("job-"):
            return {"parent_job": execution}
        elif execution.startswith("analysis-"):
            return {"parent_analysis": execution}
        else:
            raise ValueError(
                f"Search by name results in a heavy load on API server and is not currently supported by this package. "
                f"Provided execution name {execution}"
            )
