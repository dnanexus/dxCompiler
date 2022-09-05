import pytest
from concurrent.futures import ThreadPoolExecutor, as_completed
from dxcint.RegisteredTest import RegisteredTest
from copy import copy
import time
from dxcint.Messenger import State


def track_history(messenger):
    timeout = 3600
    t = 0
    history = []
    while t < timeout:
        current_state = messenger.query_state()
        history.append(current_state)
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
    results = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        h1 = executor.submit(track_history, messenger1)
        h2 = executor.submit(track_history, messenger2)
        test_futures = {
            executor.submit(test1.get_test_result),
            executor.submit(test2.get_test_result)
        }
        results = [f.result() for f in as_completed(test_futures)]
        history1 = h1.result()
        history2 = h2.result()

    assert results == [False, False]
    assert State.NOT_DONE in history1
    assert State.FINISHED in history1
    assert State.NOT_DONE in history2
    assert State.FINISHED in history2
