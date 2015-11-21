#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv -m Python/2.7.9-05st/bin -- python $0 "$@"
#}

import sys, os
import glob,shutil
import optparse
import decimal,json
from collections import OrderedDict 

def decimalEncode(obj):
    if isinstance(obj, decimal.Decimal):
        return str(obj)
    raise TypeError

# Assumes that everything after ']' or '],' and between '#' and end-of-line is a comment that belongs to the previous parameter
# obviously if this isnt the case, this function would need to be modified..
def fixupOutOfPlaceCommentsAndBadJson(filestr):
  import re
  def moveComment(matchobj):
    comment='"'+matchobj.group(2)[1:].replace('"', '').strip()+'"'
    return ', '+comment+matchobj.group(1)+'\n'

  # pull end-of-parameter-line comments into parameters
  tmp=re.sub('(\][ \t]*,?)\s+(#[^\n]*)\n', moveComment, filestr)

  # Strip full line comments
  tmp=re.sub('^[ \t]*#.*\n','',tmp,flags=re.MULTILINE)

  # Cleanup common mistake of comma after last parameter entry
  tmp=re.sub('\][ \t]*,(\s*])',r']\1',tmp)
  return tmp


def read_ff_file(file):
  fh=open(file,"r")
  filestr=fh.read()
  fh.close()
  fixStr=fixupOutOfPlaceCommentsAndBadJson(filestr)
  return json.loads(fixStr,parse_float=decimal.Decimal,object_pairs_hook=OrderedDict )

def convert_opls_to_trig(data):
  assert len(data)==5
  half=decimal.Decimal("0.5")
  zero=decimal.Decimal("0.0")
  phi,c1,c2,c3,c4 = data[:5]

  fc0 =  half*( c1+c2+c3+c4 )
  fc1 =  half*c1
  fc2 = -half*c2
  fc3 =  half*c3
  fc4 = -half*c4
  fc5 =  zero
  fc6 =  zero

  return  [phi, fc0, fc1, fc2, fc3, fc4, fc5, fc6]

def convert_and_write_params(dtype, fold, fnew):
  global toprocess
  assert dtype in toprocess

  json.encoder.ESCAPE_DCT['/']='/'
  
  pold=read_ff_file(fold)
  
  natypes,pnames=toprocess[dtype][1:]

  pnew=[]
  for row in pold:
    types=row[:natypes]
    if isinstance(row[-1],basestring):
      params=row[natypes:-1]
      memo=row[-1]
    else:
      params=row[natypes:]
      memo=''
    if("Opls" in dtype): params=convert_opls_to_trig(params)
    if(dtype == "virtuals.fda3"): params.append(decimal.Decimal("0.0"))
    #print dtype
    #print types
    #print "pnames", pnames
    #print "params", params
    if len(pnames)!=len(params):
        print "WARNING: extra elements in", dtype
        print "         Using only the first %d param columns" % len(pnames)
        print params
        print "        ", ", ".join("%s = %s" % x for x in zip(pnames,params))


    #assert len(pnames)==len(params)
    params=OrderedDict([(n,v) for n,v in zip(pnames,params)])
    rownew=OrderedDict([("type", types),("params",params),("memo",memo)])
    tmp =json.dumps(rownew,indent=None, separators=(', ', ': '),default=decimalEncode)
    pnew.append(tmp)

  s="[\n   "+",\n   ".join(pnew)+"\n]"
  fh=open(fnew,"w")
  print >>fh, s
  fh.close()

def convert_and_write_rules(fold,fnew):
  json.encoder.ESCAPE_DCT['/']='/'
  skip=["charges","pairs_es","pairs_es_balanced","pairs_lj","vdw1_14","sparsify_exclusions_and_pairs","torsion_torsion","atoms"]

  rules=read_ff_file(fold)
  pluglist=[]
  for p in rules["plugins"]:
    pname=p[0]
    pcfg=p[1]
    if pname in skip: continue
    if pname == "vdw1":
      assert pcfg == 0 or pcfg==1
    elif pname == "mass" and pcfg !=0:
      pname=p
    else:
      assert pcfg is None or pcfg==0
    pluglist.append(pname)
  rules["plugins"]=pluglist
  all=[]
  for k in rules:
    tmp=json.dumps(rules[k],indent=None, separators=(', ', ': '), default=decimalEncode)
    tmp=tmp.replace('[[','[\n'+' '*3+'[').replace(']]',']\n]').replace('],','], \n  ')
    pre='"'+k+'" : '
    tmp=pre+tmp.replace('\n','\n'+' '*(len(pre)+1))
    all.append(tmp)
  s='{\n   '+",\n".join(all).replace('\n','\n   ')+'\n}'
  fh=open(fnew,"w")
  print >>fh, s
  fh.close()      

toprocess={
  "mass"                      : ("mass",             1, ("amu",) ),
  "bonds.Harm"                : ("stretch_harm" ,    2, ("r0", "fc") ),
  "angles.Harm"               : ("angle_harm",       3, ("theta0", "fc") ),
  "ureybradley.UB"            : ("ureybradley_harm", 3, ("r0", "fc") ),
  "angles.UB"                 : ("ureybradley_harm", 3, ("r0", "fc") ),
  "propers.Proper_Trig"       : ("dihedral_trig",    4, ("phi0", "fc0", "fc1", "fc2","fc3", "fc4", "fc5", "fc6") ),
  "propers.Opls_proper"       : ("dihedral_trig",    4, ("phi0", "fc0", "fc1", "fc2","fc3", "fc4", "fc5", "fc6") ),
  "impropers.Improper_Trig"   : ("improper_trig",    4, ("phi0", "fc0", "fc1", "fc2","fc3", "fc4", "fc5", "fc6") ),
  "impropers.Opls_improper"   : ("improper_trig",    4, ("phi0", "fc0", "fc1", "fc2","fc3", "fc4", "fc5", "fc6") ),
  "impropers.Improper_Harm"   : ("improper_harm",    4, ("phi0", "fc") ),
  "vdw1.LJ12_6_sig_epsilon"   : ("vdw1",             1, ("sigma","epsilon") ),
  "vdw1_14.LJ12_6_sig_epsilon": ("vdw1_14",          1, ("sigma","epsilon") ),
  
  "torsion_torsion.cmap"      : ("torsiontorsion_cmap",8, ("cmapid",) ),
    
  "charges.BCI"               : ("charges_bci",      2, ("bci",)),
  "charges.Formal"            : ("charges_formal",   1, ("charge",) ),
  "screening.charges"         : ("charges_screening",1, ("scale12","scale13","scale14") ),
  "vdw1.exp_6x"               : ("vdw1",             1, ("alpha","epsilon","rmin")),
  "inplanewags.Harm"          : ("inplanewag_harm",  4, ("w0", "fc") ),
  "dihedrals6atom.Trig"       : ("dihedral6_trig",   4, ("phi0", "fc0", "fc2", "fc4")),
  "impropers.Improper_Anharm" : ("improper_anharm",  4, ("fc2", "fc4")),

  "virtuals.lc2"              : ("virtuals_lc2",     3, ("c1",) ),
  "virtuals.lc3"              : ("virtuals_lc3",     4, ("c1","c2") ),
  "virtuals.lc4"              : ("virtuals_lc4",     5, ("c1","c2","c3") ),    
  "virtuals.out3"             : ("virtuals_out3",    4, ("c1","c2","c3") ),
  "virtuals.fda3"             : ("virtuals_fdat3",   4, ("c1","c2","c3") ),   
  "virtuals.fdat3"            : ("virtuals_fdat3",   4, ("c1","c2","c3") ),
  "virtuals.lc2n"             : ("virtuals_lc2n",    3, ("c1",) ),
  "virtuals.lc3n"             : ("virtuals_lc3n",    4, ("c1","c2") ),
  "virtuals.lc4n"             : ("virtuals_lc4n",    5, ("c1","c2","c3") ),
  "virtuals.out3n"            : ("virtuals_out3n",   4, ("c1","c2","c3") ),
  } 
 
if __name__ == '__main__':
  
  opt = optparse.OptionParser()
  opts, args = opt.parse_args()
  if len(args) != 2:
    opt.error('Please provide an input and output file')

  dirold=args[0]
  dirnew=args[1]
  assert dirold != dirnew
  assert os.path.exists(dirold)
  assert not os.path.exists(dirnew)

  os.mkdir(dirnew)

  files=glob.glob(os.path.join(dirold,"*"))
  for file in files:
    base=os.path.basename(file)
    print file
    if base=="rules":
      convert_and_write_rules(file,os.path.join(dirnew,"rules"))
    elif base=="cmap" or base.startswith("templates") or 'read' in base.lower():
      shutil.copy(file,dirnew)
    elif base not in toprocess:
      print "I dont know how to process file: ",file
      print "File needs to be one of:\n","\n".join(sorted(toprocess.keys()))
      raise UserWarning("Invalid file: " +file)
    else:
      convert_and_write_params(base,file,os.path.join(dirnew,toprocess[base][0]))

 
