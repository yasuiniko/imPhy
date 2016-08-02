from functools import reduce
import itertools
import multiprocessing
import time

def timeit(f, s):
    bigs = s[0].upper() + s[1:]
    smalls = s[0].lower() + s[1:]
    print("{}...".format(bigs))
    t = time.time()
    x = f()
    print("Done {} in {}s.\n".format(smalls, time.time() - t))
    return x

def flatten(*lsts):
    return itertools.chain(*lsts)

def iterflatten(iterable):
    return itertools.chain.from_iterable(iterable)

def parmap(f, X, nprocs=multiprocessing.cpu_count()-1):
    """
    Taken from http://stackoverflow.com/revisions/16071616/9
    """
    def fun(f, q_in, q_out):
        while True:
            i, x = q_in.get()
            if i is None:
                break
            q_out.put((i, f(x)))

    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
            for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]

gene =           'd{}_g{}_i{}_n{}_s{}_e{}.nex'
batch_analyze =  "d{0}_g{2}_i{3}_n{5}_s{7}"
batch_general =  "d{}_g{}_i{}_n{}_s{}"
fileroot =       "d{0}_g{2}_i{3}_n{5}_s{7}_e{1}_m{4}_p{6}"
