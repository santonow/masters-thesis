from numba import njit
from multiprocessing import Process
from time import sleep


@njit
def fun(x, y):
    return x + y


def process(i, xs, ys):
    for x, y in zip(xs, ys):
        print(f"Process [{i}]", x, y, fun(x, y))
        sleep(0.2)


if __name__ == "__main__":
    procs = []
    for i in range(5):
        procs.append(
            Process(target=process, args=(i, list(range(10)), list(range(10))))
        )
        procs[-1].start()

    for proc in procs:
        proc.join()
