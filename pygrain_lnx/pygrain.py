#does stuff for coarse grained molecular dynamics, uses martini and gromacs

#import stuff
import os
import sys

class Cg():
  def __init__(self):
    '''initializes stuffs '''
    self.pydir = self.get_pydir()
    self.infiles = self.pydir + '/infiles'
    self.models = self.pydir + '/models'

  def get_pydir(self):
    '''gets the program's directory: pydir '''
    pydir = os.getcwd()
    return pydir

  def get_gro(self):
    '''gets gro files ''' 
    oldpath = '/Users/rahmadakbar/Dropbox/uds/boku2015'
    gropaths = []
    for root, dirs, files in os.walk(oldpath):
      for file in files:
        if 'prot.gro' in file:
          path = os.path.join(root, file)
          gropaths.append(path)
    return gropaths

  def cpgro(self):
    '''copies gro files to models directory '''
    gropaths = self.get_gro()
    for gropath in gropaths:
      parts = gropath.split('/')
      newparts = ['']
      for part in parts:
        if 'model' in part or 'prot.gro' in part:
          newparts.append(part)
      path = '/'.join(newparts)
      filepath = self.models + path
      dirpath = '/'.join(filepath.split('/')[:-1])
      os.system('mkdir -p %s ' % dirpath)
      os.system('cp %s %s' % (gropath, filepath))

  def get_models_gro(self):
    '''gets gro files from model directories '''
    gros = []
    for root, dirs, files in os.walk(self.pydir):
      for file in files:
        if 'prot.gro' in file:
          path = os.path.join(root, file)
          gros.append(path)
    return gros


  def findfiles(self, filename):
    filepaths = []
    for root, dirs, files in os.walk(self.pydir):
      for file in files:
        if filename in file:
          path = os.path.join(root, file)
          filepaths.append(path)
    return filepaths

  def path2names(self, path):
    '''returns path, name, namex (name with extension) from a given path '''
    parts = path.split('/')
    pathonly = '/'.join(parts[:-1])
    nameparts = parts[-1].split('.')
    name = nameparts[0]
    namex = parts[-1] # name with extension
    return pathonly, namex, name


  def do_editconf(self):
    '''converts  prot.gro files to .pdb files . Martinize takes only .pdb  '''
    prots = self.findfiles('prot.gro')
    for prot in prots[:]:
      protpath, protnamex, protname = self.path2names(prot)
      farg = protnamex
      oarg = protname + '.pdb'
      command = 'gmx editconf -f %s -o %s' % (farg, oarg)
      os.chdir(protpath)
      os.system(command)


  def do_editconf_file(self, filename):
    '''converts  prot.gro files to .pdb files . Martinize takes only .pdb  '''
    prots = self.findfiles(filename)
    for prot in prots[:]:
      filepath, namex, name = self.path2names(prot)
      farg = namex
      oarg = name + '.pdb'
      command = 'gmx editconf -f %s -o %s' % (farg, oarg)
      os.chdir(filepath)
      os.system(command)

  def do_martinize(self):
    '''executes martinize.py '''
    mpath = self.findfiles('martinize.py')[0] # martinize path
    ppaths = self.findfiles('prot.pdb') # protein paths
    for ppath in ppaths[1:]:
      filepath, namex, name = self.path2names(ppath) 
      os.chdir(filepath)
      pyarg = mpath #args are arguments
      farg = ppath
      oarg = filepath + '/%s.top' % name
      xarg = filepath + '/%s_cg.pdb' % name
      parg = 'backbone'
      ffarg = 'martini21'
      command = 'python %s -f %s -o %s -x %s -p %s -ff %s' % \
                (pyarg, farg, oarg, xarg, parg, ffarg)

      os.system(command)


  def do_martinize2(self, ppath):
    '''executes martinize.py on a protein file given in ppath '''
    mpath = self.findfiles('martinize.py')[0] # martinize path
    filepath, namex, name = self.path2names(ppath) 
    os.chdir(filepath)
    pyarg = mpath #args are arguments
    farg = ppath
    oarg = filepath + '/%s.top' % name
    xarg = filepath + '/%s_cg.pdb' % name
    parg = 'backbone'
    ffarg = 'martini21'
    command = 'python %s -f %s -o %s -x %s -p %s -ff %s' % \
              (pyarg, farg, oarg, xarg, parg, ffarg)

    os.system(command)


  def do_martinize3(self, ppath):
    '''executes martinize.py on a protein file given in ppath '''
    mpath = self.findfiles('martinize.py')[0] # martinize path
    filepath, namex, name = self.path2names(ppath) 
    os.chdir(filepath)
    pyarg = mpath #args are arguments
    farg = ppath
    oarg = filepath + '/%s.top' % name
    xarg = filepath + '/%s_cg.pdb' % name
    parg = 'backbone'
    ffarg = 'martini21'
    dssparg = self.infiles + '/dssp'
    command = 'python %s -f %s -o %s -x %s -p %s -ff %s -dssp %s' % \
              (pyarg, farg, oarg, xarg, parg, ffarg, dssparg)
    print command
    os.system(command)

  def do_insane(self):
    '''executes insane.py '''
    ipath = self.findfiles('insane.py')[0]
    ppaths = self.findfiles('prot_cg.pdb')
    for ppath in ppaths:
      filepath, namex, name = self.path2names(ppath)
      pyarg = ipath
      farg = ppath
      oarg = filepath + '/system.gro'
      parg = filepath + '/system.top'
      pbcarg = 'square'
      boxarg = '10,10,10'
      larg = 'DPPC'
      solarg = 'W'
      command = 'python %s -f %s -o %s -p %s -pbc %s -box %s -l %s -center -sol %s' \
                % (pyarg, farg, oarg, parg, pbcarg, boxarg, larg, solarg)
      os.system(command)


  def do_insane2(self, pcgpath):
    '''executes insane.py on a protein cg file given in pcgpath '''
    ipath = self.findfiles('insane.py')[0]
    filepath, namex, name = self.path2names(pcgpath)
    pyarg = ipath
    farg = pcgpath
    oarg = filepath + '/system.gro'
    parg = filepath + '/system.top'
    pbcarg = 'square'
    boxarg = '10,10,10'
    larg = 'DPPC'
    solarg = 'W'
    command = 'python %s -f %s -o %s -p %s -pbc %s -box %s -l %s -center -sol %s' \
              % (pyarg, farg, oarg, parg, pbcarg, boxarg, larg, solarg)
    os.system(command)


  def findfiles_infiles(self, filename):
    '''finds files in pydr/infiles '''
    infilepath = self.pydir + '/infiles'
    filepaths = []
    for root, dirs, files in os.walk(infilepath):
      for file in files:
        if filename in file:
          path = os.path.join(root, file)
          filepaths.append(path)
    return filepaths
        

  def include_itps(self):
    '''writes itp files onto system.top. NOTE itp orders matters, martini.itp
    must come first '''
    tops = self.findfiles('system.top')
    itps = self.findfiles_infiles('.itp')
    for top in tops:
      filepath, namex, name = self.path2names(top)
      contents = open(top).readlines()
      newcontent = ''
      for itp in itps:
        newcontent += '#include ' + '"%s"' % itp + '\n'
      protein_itp = filepath + '/Protein.itp'
      newcontent += '#include "%s"\n' % protein_itp
      for content in contents[1:]:
        newcontent += content
      newtop = open(top, 'w')
      newtop.write(newcontent)
      newtop.close()


  def include_itps2(self,toppath):
    '''writes itp files onto a topology file ( system.top). NOTE itp orders matters, martini.itp
    must come first '''
    itps = self.findfiles_infiles('.itp')
    itps = sorted(itps)
    filepath, namex, name = self.path2names(toppath)
    contents = open(toppath).readlines()
    newcontent = ''
    for itp in itps:
      newcontent += '#include ' + '"%s"' % itp + '\n'
    for content in contents[1:]:
      newcontent += content
    newtop = open(toppath, 'w')
    newtop.write(newcontent)
    newtop.close()


  def grompp(self, farg, parg, carg, oarg):
    ''' executes gmx grompp'''
    command = 'gmx grompp -f %s -p %s -c %s -o %s ' \
              % ( farg, parg, carg, oarg)
    os.system(command)
    print '\n\n' + command + '\n\n'
    print os.getcwd()

  def grompp_ndx(self, farg, parg, carg, oarg, narg):
    ''' executes gmx grompp'''
    command = 'gmx grompp -f %s -p %s -c %s -o %s -n %s ' \
              % ( farg, parg, carg, oarg, narg)
    os.system(command)
    print '\n\n' + command + '\n\n'
    print os.getcwd()

  def nt1mdrun(self, tprname):
    '''executes gmx mdrun number of thread is 1 '''
    command = 'gmx mdrun -deffnm %s -nt 1 -v' % tprname
    os.system(command)
    print '\n\n' + command + '\n\n'
    print os.getcwd()

  def do_genion(self, ionsmdp, systemfile):
    '''adds ions to a systemfile '''
    systempath, systemnamex, systemname = self.path2names(systemfile)
    farg = ionsmdp
    carg = systemnamex
    parg = systemname + '.top'
    oarg = 'ions.tpr'
    self.grompp(farg, parg, carg, oarg)
    sarg = oarg
    oarg2 = systemnamex
    pnamearg = 'NA' 
    nnamearg = 'CL'
    command = 'echo W | gmx genion -s %s -o %s -p %s -pname %s -nname %s -neutral -conc 0.1' \
              % (sarg, oarg2, parg, pnamearg, nnamearg)
    os.system(command) 

  def do_minim(self, minimmdp, systemfile):
    '''does minimization '''
    systempath, systemnamex, systemname = self.path2names(systemfile)
    farg = minimmdp
    parg = systemname + '.top'
    carg = systemnamex
    oarg = 'minim.tpr'
    self.grompp(farg, parg, carg, oarg)
    self.nt1mdrun('minim')


  def do_minim_mult(self, minimmdp, systemfile, n):
    '''does minimization '''
    systempath, systemnamex, systemname = self.path2names(systemfile)
    farg = minimmdp
    parg = systemname + '.top'
    carg = systemnamex
    oarg = 'minim.tpr'
    self.grompp(farg, parg, carg, oarg)
    self.nt1mdrun('minim')
    for i in range(n):
      farg = minimmdp
      parg = systemname + '.top'
      carg = 'minim.gro'
      oarg = 'minim.tpr'
      self.grompp(farg, parg, carg, oarg)
      self.nt1mdrun('minim')
   

  def make_ndx(self, farg, oarg):
    '''executes gmx make_ndx '''
    command = 'echo "1|13 \n 14|15 \n q" | gmx make_ndx -f %s -o %s' % (farg, oarg)    
    os.system(command)

  def do_nvt(self, nvtmdp, systemfile):
    '''does minimization '''
    systempath, systemnamex, systemname = self.path2names(systemfile)
    self.make_ndx('minim.gro', 'nvt.ndx') 
    farg = nvtmdp
    parg = systemname + '.top'
    carg = 'minim.gro'
    oarg = 'nvt.tpr'
    narg = 'nvt.ndx'
    self.grompp_ndx(farg, parg, carg, oarg, narg)
    self.nt1mdrun('nvt')


  def do_npt(self, nptmdp, systemfile):
    '''does minimization '''
    systempath, systemnamex, systemname = self.path2names(systemfile)
    #self.make_ndx('minim.gro', 'nvt.ndx') 
    farg = nptmdp
    parg = systemname + '.top'
    carg = 'nvt.gro'
    oarg = 'npt.tpr'
    narg = 'nvt.ndx'
    self.grompp_ndx(farg, parg, carg, oarg, narg)
    self.nt1mdrun('npt')


  def prepare_martini(self):
    '''prepares starting pdb files. uses martinize and insane. '''
    self.do_martinize()
    self.do_insane()
    self.include_itps()  

  def prepare_gmx(self):
    '''gromacs. prepares simulation, uses the given mdp files '''
    minimmdp  = self.findfiles('minim.mdp')[0]
    nvtmdp = self.findfiles('nvt.mdp')[0]
    nptmdp = self.findfiles('npt.mdp')[0]
    systems = self.findfiles('system.gro')
    for system in systems[1:]:
      systempath, systemnamex, systemname = self.path2names(system)
      os.chdir(systempath)
      self.do_minim(minimmdp, system)
      self.do_nvt(nvtmdp, system)
      self.do_npt(nptmdp, system)

  def reset_models(self):
    '''reset model directories, leaves only prot.gro and prot.pdb '''
    modeldir = self.pydir + '/models'
    wanteds = ['prot.gro', 'prot.pdb']
    for root, dirs, files in os.walk(modeldir):
      for file in files:
        if file not in wanteds:
          filepath = os.path.join(root, file)
          os.system('rm %s' % filepath)

  def prepare_martini_gmx(self):
    '''takes prot.pdb files, executes the following: martini, insane, minim, nvt,
    and npt '''
    protfiles  = self.findfiles('prot.pdb') 
    for protfile in protfiles[:]:
      protpath, protnamex, protname = self.path2names(protfile)
      dircontents1 = os.listdir(protpath)
      print dircontents1
      nptfound = False
      if 'npt.gro' not in dircontents1:
        nptfound = False
      else:
        nptfound = True
      while nptfound == False:
        print protpath
        self.do_martinize2(protfile)
        pcgpath = protpath + '/prot_cg.pdb'
        self.do_insane2(pcgpath)
        toppath = protpath + '/system.top'     
        self.include_itps2(toppath)
        minimmdp  = self.findfiles('minim.mdp')[0]
        nvtmdp = self.findfiles('nvt.mdp')[0]
        nptmdp = self.findfiles('npt.mdp')[0]
        systemfile = protpath + '/system.gro'
        os.chdir(protpath)
        self.do_minim(minimmdp, systemfile)
        self.do_nvt(nvtmdp, systemfile)
        self.do_npt(nptmdp, systemfile)
        dircontents2 = os.listdir(protpath)
        if 'npt.gro' in dircontents2:
          nptfound = True
          print 'npt done for %s' % protfile

  def merge_martini_insane_top(self, martinitop, insanetop):
    '''merges martini and insane topology files'''
    mtop = open(martinitop).readlines()
    itop = open(insanetop).readlines()
    newtop = mtop + ['\n'] +  itop[-3:]
    newcontent = ''.join(newtop)
    newfile = open(insanetop, 'w')
    newfile.write(newcontent)
    newfile.close()

  def prepare_martini_gmx_multminim(self, n):
    '''takes prot.pdb files, executes the following: martini, insane, minim, nvt,
    and npt '''
    protfiles  = self.findfiles('prot.pdb') 
    for protfile in protfiles[:]:
      protpath, protnamex, protname = self.path2names(protfile)
      dircontents1 = os.listdir(protpath)
      nptfound = False
      if 'npt.gro' not in dircontents1:
        nptfound = False
      else:
        nptfound = True
      while nptfound == False:
        print protpath
        self.do_martinize3(protfile)
        pcgpath = protpath + '/prot_cg.pdb'
        self.do_insane2(pcgpath)
        toppath = protpath + '/system.top'
        martinitop = protpath + '/prot.top' 
        insanetop = protpath + '/system.top' 
        self.merge_martini_insane_top(martinitop, insanetop)  
        self.include_itps2(toppath)
        print protfile
        sys.exit()
        minimmdp  = self.findfiles('minim.mdp')[0]
        nvtmdp = self.findfiles('nvt.mdp')[0]
        nptmdp = self.findfiles('npt.mdp')[0]
        systemfile = protpath + '/system.gro'
        os.chdir(protpath)
        self.do_minim_mult(minimmdp, systemfile, n)
        #sys.exit()
        self.do_nvt(nvtmdp, systemfile)
        self.do_npt(nptmdp, systemfile)
        dircontents2 = os.listdir(protpath)
        if 'npt.gro' in dircontents2:
          nptfound = True
          print 'npt done for %s' % protfile

  def prepare_martini_insane(self):
    '''prepares martini and insane stuff '''
    protfiles  = self.findfiles('prot.pdb') 
    for protfile in protfiles[:]:
      protpath, protnamex, protname = self.path2names(protfile)
      self.do_martinize3(protfile)
      pcgpath = protpath + '/prot_cg.pdb'
      self.do_insane2(pcgpath)
      toppath = protpath + '/system.top'
      martinitop = protpath + '/prot.top' 
      insanetop = protpath + '/system.top' 
      self.merge_martini_insane_top(martinitop, insanetop)  
      self.include_itps2(toppath)
      print protfile

  def gmx_minim_eq(self):
    '''executes the following minim, nvt,
    and npt '''
    protfiles  = self.findfiles('prot.pdb') 
    for protfile in protfiles[:]:
      protpath, protnamex, protname = self.path2names(protfile)
      dircontents1 = os.listdir(protpath)
      nptfound = False
      if 'npt.gro' not in dircontents1:
        nptfound = False
      else:
        nptfound = True
      while nptfound == False:
        minimmdp  = self.findfiles('minim.mdp')[0]
        nvtmdp = self.findfiles('nvt.mdp')[0]
        nptmdp = self.findfiles('npt.mdp')[0]
        systemfile = protpath + '/system.gro'
        os.chdir(protpath)
        self.do_minim(minimmdp, systemfile)
        sys.exit()
        self.do_nvt(nvtmdp, systemfile)
        self.do_npt(nptmdp, systemfile)
        dircontents2 = os.listdir(protpath)
        if 'npt.gro' in dircontents2:
          nptfound = True
          print 'npt done for %s' % protfile



  def prepare_martini_gmx2(self):
    '''takes prot.pdb files, executes the following: martini, insane, genion, minim, nvt,
    and npt '''
    protfiles  = self.findfiles('prot.pdb') 
    for protfile in protfiles[:]:
      protpath, protnamex, protname = self.path2names(protfile)
      dircontents1 = os.listdir(protpath)
      print dircontents1
      nptfound = False
      if 'npt.gro' not in dircontents1:
        nptfound = False
      else:
        nptfound = True
      while nptfound == False:
        print 'not yet'
        print protpath
        self.do_martinize2(protfile)
        pcgpath = protpath + '/prot_cg.pdb'
        self.do_insane2(pcgpath)
        toppath = protpath + '/system.top'     
        self.include_itps2(toppath)
        ionsmdp  = self.findfiles('ions.mdp')[0]
        minimmdp  = self.findfiles('minim.mdp')[0]
        nvtmdp = self.findfiles('nvt.mdp')[0]
        nptmdp = self.findfiles('npt.mdp')[0]
        systemfile = protpath + '/system.gro'
        os.chdir(protpath)
        self.do_genion(ionsmdp, systemfile)
        self.do_minim(minimmdp, systemfile)
        self.do_nvt(nvtmdp, systemfile)
        #sys.exit('out')
        self.do_npt(nptmdp, systemfile)
        print 'npt done for %s' % protfile
        dircontents2 = os.listdir(protpath)
        if 'npt.gro' in dircontents2:
          nptfound = True
          print 'npt done for %s' % protfile



  def add_chain(self):
    '''adds chain information to prot.pdb files '''
    pdbfiles = self.findfiles('prot.pdb')
    for pdbfile in pdbfiles:
      contents = open(pdbfile).readlines()
      newcontents = ''
      for content in contents:
        if 'ATOM' not in content[:5]:
          newcontents += content
        else:
          serial = int(content[6:11])
          if serial < 270:
            newcontent = content[:21] + 'S' + content[22:]
            newcontents += newcontent
            counter = serial
          elif serial == 270:
            newcontent = content[:21] + 'S' + content[22:]
            newcontents += newcontent
            newcontents += 'TER\n'
          else:
            newcontent = content[:21] + 'A' + content[22:]
            newcontents += newcontent
      newpdb = open(pdbfile, 'w')
      newpdb.write(newcontents)
      newpdb.close()


  def check_system(self):
    '''returns basics plots '''
    edrs = self.findfiles('.edr')
    edrs = [item for item in edrs if '#' not in item]
    metrics = ['Potential', 'Pressure', 'Temperature']
    for edr in edrs:
      edrdir, edrnamex, edrname = self.path2names(edr)
      farg = edr
      for metric in metrics:
        oarg = edrdir + '/' + edrname + '_' + metric
        command = 'echo "%s \n" | gmx energy -f %s -o %s' % (metric, farg, oarg)
        os.system(command)
  


  def translate_x(self, addx):
    '''translates x cordinates by adding additonal x(addx). '''
    pdbfiles = self.findfiles('prot.pdb')
    for pdbfile in pdbfiles:
      contents = open(pdbfile).readlines()
      newcontents = ''
      for content in contents:
        if 'ATOM' not in content[:5]:
          newcontents += content
        else:
          serial = int(content[6:11])
          if serial <= 270:
            xcord = content[30:38]
            newxcord = float(xcord) + addx
            newxcord = '%8s' % newxcord
            newcontent = content[:30] + newxcord + content[38:]
            newcontents += newcontent
          else:
            newcontents += content
      pdbdir, namex, name = self.path2names(pdbfile)
      pdbold = pdbdir + '/' + name + '.old'
      os.system('cp %s %s' % (pdbfile, pdbold))
      newpdb = open(pdbfile, 'w')
      newpdb.write(newcontents)
      newpdb.close()


### cg stuff ###
  def cg_wf(self):
    '''workflow for cgmd '''
    #self.do_editconf()
    #self.add_chain()
    #self.do_martinize()
    #self.do_insane()
    #self.include_itps()
    #self.prepare_martini()
    #self.prepare_gmx()
    #self.reset_models()
    #self.prepare_martini_gmx()
    #self.prepare_martini_gmx_multminim(0)
    #self.check_system()
    #self.translate_x(50)
    #self.do_editconf_file('system.gro')
    #self.prepare_martini_insane()
    #self.gmx_minim_eq()

### end cg stuff ###

###dev notes###
#c = Cg()
#c.cg_wf()
