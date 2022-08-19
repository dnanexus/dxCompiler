import logging
import dxpy

from typing import Optional, Dict


class MessengerError(Exception):
    """
    Class to handle Messenger Errors
    """


class Status(object):
    NOT_STARTED, RUNNING, FINISHED, FAIL = "NOT_STARTED", "RUNNING", "FINISHED", "FAIL"
    _allowed_states = {NOT_STARTED, RUNNING, FINISHED, FAIL}

    def __init__(self):
        """
        Class to maintain the status of current test
        """
        self._state = self.NOT_STARTED

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, new_state: str):
        if new_state not in self._allowed_states:
            raise MessengerError(
                f"Status.state(): state {new_state} is not allowed for this messenger. Allowed states are "
                f"{self._allowed_states}"
            )
        self._state = new_state


class Messenger(object):
    def __init__(self, test_name: str, job_id: str, variant: Optional[int] = None):
        self._test_name = test_name
        self._variant = variant if variant else ""
        self._status = Status()
        self._status.state = Status.RUNNING
        self._job = dxpy.DXJob(job_id)
        self._interval = 20     # Default job polling interval

    @property
    def state(self):
        return self._status.state

    def wait_for_completion(self) -> str:
        logging.info(f"Waiting for completion of the test `{self._test_name}` with input variant {self._variant}")
        test_name_with_variant = ".".join([self._test_name, self._variant]).strip(".")
        try:
            self._job.wait_on_done(interval=self._interval)
            logging.info(f"Analysis {test_name_with_variant} succeeded")
            self._status.state = Status.FINISHED
        except dxpy.DXJobFailureError:
            logging.info(f"Analysis {test_name_with_variant} failed. This failure will verified")
            self._status.state = Status.FAIL
        return self.state
