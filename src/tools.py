import time

def timeit(f, s):
    bigs = s[0].upper() + s[1:]
    smalls = s[0].lower() + s[1:]
    print("{}...".format(bigs))
    t = time.time()
    x = f()
    print("Done {} in {}s.\n".format(smalls, time.time() - t))
    return x

gene =     'd{}_g{}_i{}_n{}_s{}_e{}.nex'
batch =    "d{0}_g{2}_i{3}_n{5}_s{7}"
fileroot = "d{0}_g{2}_i{3}_n{5}_s{7}_e{1}_m{4}_p{6}"
