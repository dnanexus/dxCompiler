import logging
import dxpy

from typing import Optional, Union
from enum import Enum, auto


class State(Enum):
    NOT_DONE = auto()  # 'idle' or "runnable" or "running" or "in_progress"
    FINISHED = auto()  # 'done'
    FAIL = auto()  # 'failed' or "terminated" or "partially_failed", "terminating"
    UNKNOWN = (
        auto()
    )  # 'unknown' or "waiting_on_input", "waiting_on_output" "debug_hold"


class Messenger(object):
    """Class for working with the state of a analysis/job."""

    def __init__(
        self,
        test_name: str,
        job_id: str,
        variant: Optional[object] = None,
        interval: int = 20,
    ):
        self._test_name = test_name
        self._variant = variant if variant else ""
        self._state: State = State.NOT_DONE
        self._execution: Union[
            dxpy.DXAnalysis, dxpy.DXJob
        ] = dxpy.bindings.dxdataobject_functions.get_handler(job_id)
        self._interval = interval

    @property
    def state(self) -> State:
        return self._state

    def wait_for_completion(self) -> State:
        """Waits for the job to finish and returns its state.

        It uses a blocking call to dxpy, but can be used concurrently in a separate thread.
        Most common usage is in the _validate method of RegisteredTest descendant.

        Returns:
            State.FINISHED if the execution ran to completion, State.FAIL if it failed or was terminated.
        """
        logging.info(
            f"Waiting for completion of the test `{self._test_name}` with input variant {self._variant}"
        )
        test_name_with_variant = f"{self._test_name}.{self._variant}"
        try:
            self._execution.wait_on_done(interval=self._interval)
            logging.info(f"Analysis {test_name_with_variant} succeeded")
            self._state = State.FINISHED
        except dxpy.DXJobFailureError:
            logging.info(
                f"Analysis {test_name_with_variant} failed. This failure will verified"
            )
            self._state = State.FAIL
        return self.state

    def query_state(self, extended: bool = False) -> Union[State, str]:
        """Queries the state of the execution and returns it in the form of a State enum or a string.
        Args:
            extended: bool. If True, returns the state as a string from the API, otherwise a simple State enum.
        Returns:
            Union[State,str] State enum or a string state from the API.
        """
        returned_state = self._execution._get_state()
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
