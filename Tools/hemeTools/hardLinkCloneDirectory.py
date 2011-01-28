import os.path

def clone(sourceDir, destDir):
    assert os.path.exists(sourceDir), "Source directory does not exist"
    assert os.path.isdir(sourceDir), "Source is not a directory"
    assert not os.path.exists(destDir), "Destination already exists"

    os.mkdir(destDir)

    filesToLink = ['config.dat', 'config.xml', 'coords.asc',
                   'pars.asc', 'rt_pars.asc']

    for f in filesToLink:
        src = os.path.join(sourceDir, f)
        assert os.path.exists(src), 'Source file "%s" does not exist' % src
        dest = os.path.join(destDir, f)
        os.link(src, dest)
        continue
    
    return

if __name__ == '__main__':
    import sys
    
    sourceDir = sys.argv[1]
    destDir = sys.argv[2]
    clone(sourceDir, destDir)

    
