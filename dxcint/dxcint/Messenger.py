import dxpy

from typing import Optional, Union, Dict
from enum import Enum, auto

from dxcint.Context import Context


class State(Enum):
    NOT_DONE = auto()  # 'idle' or "runnable" or "running" or "in_progress"
    FINISHED = auto()  # 'done'
    FAIL = auto()  # 'failed' or "terminated" or "partially_failed", "terminating"
    UNKNOWN = (
        auto()
    )  # 'unknown' or "waiting_on_input", "waiting_on_output" "debug_hold"


class Messenger(object):
    def __init__(
        self,
        context: Context,
        test_name: str,
        job_id: str,
        variant: Optional[int] = None,
        interval: int = 20
    ):
        """
        Class to communicate between RegisteredTest and the respective analysis/job on the platform.
        Args:
            context: instance of dxcint.Context
            test_name: name of the test as it's been registered in the suite
            job_id: job or analysis id of the running test on the platform
            variant: variant of the input if multiple inputs were used for the same test workflow
            interval: Number of seconds between queries to the job’s state
        """
        self._context = context
        self._test_name = test_name
        self._variant = variant if variant else ""
        self._state: State = State.NOT_DONE
        self._execution: Union[
            dxpy.DXAnalysis, dxpy.DXJob
        ] = dxpy.bindings.dxdataobject_functions.get_handler(job_id)
        self._interval = interval
        self._execution_describe = None

    @property
    def state(self) -> State:
        return self._state

    @property
    def describe(self) -> Dict:
        if not self._execution_describe:
            self._execution_describe = self._execution.describe()
        return self._execution_describe

    def wait_for_completion(self) -> State:
        """
        Waits for the job to finish and returns its state.
        It uses a blocking call to dxpy, but can be used concurrently in a separate thread.
        Most common usage is in the _validate method of an instance of the RegisteredTest subclass.

        Returns:
            State.FINISHED if the execution ran to completion, State.FAIL if it failed or was terminated.
        """
        self._context.logger.info(
            f"Waiting for completion of the test `{self._test_name}` with input variant {self._variant}"
        )
        test_name_with_variant = f"{self._test_name}.{self._variant}"
        try:
            self._execution.wait_on_done(interval=self._interval)
            self._context.logger.info(f"Analysis {test_name_with_variant} succeeded")
            self._state = State.FINISHED
        except dxpy.DXJobFailureError:
            self._context.logger.info(
                f"Analysis {test_name_with_variant} failed. This failure will verified"
            )
            self._state = State.FAIL
        return self.state

    def query_state(self, extended: bool = False) -> Union[State, str]:
        """
        Queries the state of the execution and returns it in the form of a State enum or a string.
        Args:
            extended: bool. If True, returns the state as a string from the API, otherwise a simple State enum.
        Returns:
            Union[State,str] State enum or a string state from the API.
        """
        returned_state = self._execution.describe(fields=dict(state=True))["state"]
        if extended:
            self._state = returned_state
        elif returned_state in [
            "terminated",
            "failed",
            "partially_failed",
            "terminating",
        ]:
            self._state = State.FAIL
        elif returned_state in ["done"]:
            self._state = State.FINISHED
        elif returned_state in ["idle", "in_progress", "runnable", "running"]:
            self._state = State.NOT_DONE
        else:
            self._state = State.UNKNOWN

        return self.state

    def describe_execution(self) -> Dict:
        """Returns the description of the execution."""
        return self._execution.describe()
