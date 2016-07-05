import sys
import os

if __name__ == "__main__":
    """
    Converts Mesquite nexus file to R readable nexus files. R files are 
    written to r_nexus/. The folder is created if it does not exist.

    Usage: python3 mesquite_to_R.py fdir

    fdir is the directory holding the Mesquite nexus files or holding
    a directory called mesquite_nexus, which contains the nexus files.
    """

    # get the folder name
    if len(sys.argv) > 1:
        fdir = sys.argv[1]
    else:
        raise NameError("Please provide a folder to read.")

    # organize filesystem
    if 'mesquite_nexus' in os.listdir(fdir):
        r_dir = os.path.join(fdir, "r_nexus")
        fdir = os.path.join(fdir, 'mesquite_nexus')
    else:
        r_dir = os.path.join(os.path.dirname(fdir), "r_nexus")
    if not os.path.isdir(r_dir):
        os.makedirs(r_dir)

    # conversion loop
    for fname in os.listdir(fdir):

        # only process .nex files
        if not fname[-4:] == ".nex":
            print("Skipping file {} since it".format(fname) + 
                  " does not end in .nex.")
            continue
        
        # read file
        fpath = os.path.join(fdir, fname)
        with open(fpath, 'r') as f:
            in_tree = False
            in_associates = False
            associates = []
            fcontents = []

            for line in f:

                # identify species/individuals associations
                if in_associates:
                    assoc = line.split('/')

                    if len(assoc) > 1:
                        species = assoc[0].strip()
                        inds = assoc[1].replace(',', '').strip().split(' ')
                        associates.append([species]+inds)

                    if ';' in line:
                        in_associates = False
                        continue

                if "ASSOCIATES" in line:
                    in_associates = True

                # identify trees
                if "BEGIN TREES" in line:
                    in_tree = True
                    line = "BEGIN TREES;\n"

                if in_tree:
                    fcontents.append(line)

                    if "Title" in line:
                        tree_name = line.split("Title ")[-1][:-2]

                    if "END" in line:
                        in_tree = False

                        # write fcontents to file
                        outname = "{}_{}.nex".format(fname[:-4], tree_name)
                        outpath = os.path.join(r_dir, outname)
                        with open(outpath, 'w') as g:
                            for line in associates:
                                [g.write(x+", ") for x in line]
                                g.write('\n')
                            for line in fcontents:
                                g.write(line)

                        fcontents = []

