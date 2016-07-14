import time

def timeit(f, s):
    bigs = s[0].upper() + s[1:]
    smalls = s[0].lower() + s[1:]
    print("{}...".format(bigs))
    t = time.time()
    x = f()
    print("Done {} in {}s.\n".format(smalls, time.time() - t))
    return x