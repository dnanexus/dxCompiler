import pytest
from concurrent.futures import ThreadPoolExecutor, as_completed
from dxcint.RegisteredTest import RegisteredTest
from copy import copy
from typing import List, Tuple
import time
from dxcint.Messenger import State


def track_history(messenger) -> List[Tuple[float, State]]:
    timeout = 3600
    t = 0
    history = []
    while t < timeout:
        current_state = messenger.query_state()
        history.append((time.time(), current_state))
        if current_state == State.FINISHED or current_state == State.FAIL:
            break
        time.sleep(15)
        t += 15
    return history


@pytest.mark.usefixtures("build_executable_wdl")
def test_parallel_execution(build_executable_wdl: RegisteredTest):
    test1 = build_executable_wdl
    test2 = copy(build_executable_wdl)
    messenger1 = test1.messenger
    messenger2 = test2.messenger
    with ThreadPoolExecutor(max_workers=10) as executor:
        h1 = executor.submit(track_history, messenger1)
        h2 = executor.submit(track_history, messenger2)
        test_futures = {
            executor.submit(test1.get_test_result),
            executor.submit(test2.get_test_result),
        }
        results = [f.result() for f in as_completed(test_futures)]
        history1 = h1.result()
        history2 = h2.result()

    assert results == [False, False]
    assert State.NOT_DONE in [state for _, state in history1]
    assert State.FINISHED in [state for _, state in history1]
    assert State.NOT_DONE in [state for _, state in history2]
    assert State.FINISHED in [state for _, state in history2]

    # the histories should share a period at which both tests are not done to confirm concurrency
    # (stronger test is not possible, because the first job might finish before the second one gets a worker)
    assert any(
        abs(ts1 - ts2) < 10 and state1 == State.NOT_DONE and state2 == State.NOT_DONE
        for ts1, state1 in history1
        for ts2, state2 in history2
    )
