import cpsfile

derfile=cpsfile.derFile('/home/leon/code/CPSPy/test_modesum_ti/TRDER.TXT')
tree= derfile.get_tree()
derfile2=cpsfile.derFile('/home/leon/code/CPSPy/test_modesum_ti/TLDER.TXT')
tree2 = derfile2.get_tree(tree)