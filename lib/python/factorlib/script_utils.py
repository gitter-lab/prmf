import os, os.path
import sys
import subprocess as sp
import multiprocessing as mp
import networkx as nx
import distutils.spawn
import hashlib

def run_command(outdir, cmd, *args, **kwargs):
  """Run command and throw error if non-zero exit code

  Keyword Arguments
  -----------------
  condor : boolean
  stdout_fh : io-like
  stderr_fh : io-like

  TODO
  ----
  change interface to match format_vars
  """
  if not 'condor' in kwargs:
    kwargs['condor'] = False
  args = [cmd] + list(args)
  if(kwargs['condor']):
    args = ["condor_submitter.sh"] + args

  stdout_fh = None
  if not 'stdout' in kwargs:
    stdout_fh = open(os.path.join(outdir, "{}.out".format(cmd)), "w")
  else:
    stdout_fh = kwargs['stdout']

  stderr_fh = None
  if not 'stderr' in kwargs:
    stderr_fh = open(os.path.join(outdir, "{}.err".format(cmd)), "w")
  else:
    stderr_fh = kwargs['stderr']

  sys.stdout.write("[STATUS] Launching {}\n".format(str(args)))
  complete_proc = sp.check_call(args, stdout=stdout_fh, stderr=stderr_fh)

def run_command_cp(outdir, cmd, *args, **kwargs):
  """
  "Decorator" around run_command to enable checkpointing
  """
  if 'no_checkpoint' in kwargs:
    kwargs['no_checkpoint'] = True
  else:
    kwargs['no_checkpoint'] = False

  if kwargs['no_checkpoint']:
    run_command(outdir, cmd, *args, **kwargs)
  else:
    # identify this command and its arguments (TODO and its environment?) with a hash value
    to_hash = [cmd] + list(args)
    hash_v = hash("".join(to_hash))
    hash_fp = os.path.join(outdir, str(hash_v))

    # check if command has previously been run successfully
    if os.path.exists(hash_fp):
      # if so, dont run again
      pass
    else:
      run_command(outdir, cmd, *args, **kwargs)

      # if no error from check_call, command successful, write a file indicating this
      with open(hash_fp, 'w') as fh:
        fh.write("\t".join(to_hash))

def get_condor_submit_fp():
  """
  TODO
  ----
  allow a different environment variable to specify the location of this file
  """
  home_dir = os.environ.get("HOME")
  if(home_dir is None):
    raise ValueError("Required environment variable: HOME")
  return os.path.join(home_dir, ".condor", "submitter.sub")

def format_vars(job_id, exe=None, args=[], out=None, err=None, requirements=None, **kwargs):
  """
  Parameters
  -------
  job_id : str
    Condor DAG job identifier

  exe : str
    executable to run

  args : list of str
    arguments to exe

  out : str
    filepath to write stdout to

  err : str
    filepath to write stderr to

  requirements : str
    Condor-style requirements statement e.g. 'OpSysMajorVer == 7'

  Returns
  -------
  vars_stmt : str
  """
  VARS_FMT_STR = "VARS {} executable=\"{}\""
  exe_path = get_exe_path(exe)
  vars_stmt = VARS_FMT_STR.format(job_id, exe_path, " ".join(args))
  if(exe is None):
    raise ValueError("keyword argument \"exe\" is required")
  if(len(args) > 0):
    vars_stmt += " arguments=\"{}\"".format(" ".join(args))
  if(out is not None):
    vars_stmt += " output=\"{}\"".format(out)
  if(err is not None):
    vars_stmt += " error=\"{}\"".format(err)
  if(requirements is not None):
    vars_stmt += "requirements = \"{}\"".format(requirements)
  return vars_stmt

def get_exe_path(inpath):
  """
  Wrapper around distutils.spawn.find_executable to provide warning if found executable is in cwd
  """
  outpath = os.path.abspath(distutils.spawn.find_executable(inpath))
  cwd_exe_path = os.path.abspath(os.path.join(os.curdir, inpath))
  if(outpath == cwd_exe_path):
    # this exe may surprise some users because find_executable behaves differently from "which" in
    # this respect: find_executable will also include the current working directory in the search
    # for the executable
    sys.stderr.write("[warning] found executable {} is in current working directory and may be a different version of {} than found according to the command-line program \"which\"\n".format(outpath, exe))
  return outpath

def job_attrs_to_job_name(exe=None, args=None, out=None, err=None, **kwargs):
  """
  http://research.cs.wisc.edu/htcondor/manual/v8.7/2_10DAGMan_Applications.html#SECTION003102100000000000000
  The JobName can be any string that contains no white space, except for the strings PARENT 
  and CHILD (in upper, lower, or mixed case). JobName also cannot contain special characters 
  ('.', '+') which are reserved for system use.
  """
  hash_obj = hashlib.sha256()
  hash_obj.update(exe.encode())
  hash_obj.update(" ".join(args).encode())
  job_name = hash_obj.hexdigest()
  return job_name

def write_condor_dag(dag_fp, digraph):
  """
  See run_digraph
  """
  JOB_FMT_STR = "JOB {} {}"
  PARENT_FMT_STR = "PARENT {} CHILD {}"
  int_to_chr_offset = ord('A')
  job_int_to_name = {}

  with open(dag_fp, "w") as fh:
    # get generic condor job description filepath
    condor_submit_fp = get_condor_submit_fp()

    job_order = nx.topological_sort(digraph)
    for job_int in job_order:
      # write JOB declaration for root
      job_int = root
      job_attrs = digraph.node[job_int]
      job_name = job_attrs_to_job_name(**job_attrs)
      job_int_to_name[job_int] = job_name
      fh.write(JOB_FMT_STR.format(job_name, condor_submit_fp) + "\n")

      # write VARS declaration for root
      fh.write(format_vars(job_name, **job_attrs) + "\n")

    # write all edges at end
    for edge in digraph.edges_iter():
      # write PARENT .. CHILD declaration
      fh.write(PARENT_FMT_STR.format(job_int_to_name[edge[0]], job_int_to_name[edge[1]]) + "\n")

def submit_condor_dag(dag_fp):
  # TODO output files go in cwd or location of dag file
  args = ["condor_submit_dag", dag_fp]
  sp.check_call(args)

def bfs_nodes(G, source):
  node_list = []
  edge_list = list(nx.bfs_edges(G, source))
  node_list.append(edge_list[0][0])
  for edge in edge_list:
    node_list.append(edge[1])
  return node_list

def run_digraph(outdir, digraph, condor=False, dry_run=False, root_node=0, exit_on_err=True, **kwargs):
  """
  Run a set of jobs specified by a directed (acyclic) graph

  Parameters
  ----------
  digraph : nx.DiGraph
    directed graph specifying job execution order; nodes in the graph are assumed to have node
    attribute 'args' which specifies the executable in args[0] and its arguments in args[1:] and node identifiers
    in [0..n]

  condor : bool
    if True, use condor_submitter.sh to submit jobs

  dry_run : bool
    if True, do not submit the DAG

  exit_on_err : bool
    if True (and condor False), stop execution of digraph when one of the job nodes fails

  Returns : TODO
  -------
  job_ids : list of str
  """
  job_ids = []
  if(condor):
    # then create a DAG description file for Condor
    dag_fp = os.path.join(outdir, "digraph.dag")
    # TODO job_ids?
    write_condor_dag(dag_fp, digraph)
    if(not dry_run):
      submit_condor_dag(dag_fp)
  else:
    # TODO even with the pool it appears only 1 process is running
    pool = mp.Pool(processes=mp.cpu_count()-1)
    digraph = digraph.copy() # add proc node attr pointing to a Popen object
    #job_order = bfs_nodes(digraph, root_node)
    job_order = nx.topological_sort(digraph)
    for job_id in job_order:
      # wait for any predecessors to finish
      #preds = digraph.predecessors(job_id)
      #for pred in preds:
      #  digraph.node[pred]['async_result'].wait()

      #sum_v = 0
      #for exit in pred_exits:
      #  sum_v += exit
      #if(sum_v > 0):
      #  # then a predecessor failed
      #  raise RuntimeError("[ERROR] a predecessor to {} failed".format(digraph.node[job_id]['exe']))

      # launch this node's job
      job_attrs = digraph.node[job_id]
      args = [job_attrs['exe']] + job_attrs['args']
      stdout_fh = None
      if('out' in job_attrs):
        stdout_fh = open(job_attrs['out'], 'w')
      else:
        stdout_fh = sys.stdout
      stderr_fh = None
      if('err' in job_attrs):
        stderr_fh = open(job_attrs['err'], 'w')
      else:
        stderr_fh = sys.stderr
      sys.stdout.write("[STATUS] Launching {} > {} 2> {}\n".format(" ".join([digraph.node[job_id]['exe']] + digraph.node[job_id]['args']), job_attrs['out'], job_attrs['err']))

      #def callback_for_job(exit_status):
      #  digraph.node[job_id]['exit'] = exit_status
      #digraph.node[job_id]['exit'] = None
      # TODO for some reason providing stdout and stderr breaks everything?
      #async_result = pool.apply_async(sp.check_call, [args]) #, {'stdout': stdout_fh, 'stderr': stderr_fh})
      #digraph.node[job_id]['async_result'] = async_result
      #async_result.wait()

      # TODO synchronous only
      #proc = sp.Popen(args, stdout=stdout_fh, stderr=stderr_fh)
      exit_code = -1
      if exit_on_err:
        exit_code = sp.check_call(args, stdout=stdout_fh, stderr=stderr_fh)
      else: 
        exit_code = sp.call(args, stdout=stdout_fh, stderr=stderr_fh)
      print(exit_code)

  return job_ids

def parse_condor_submit_stdout(stdout):
  """
  Returns
  -------
  job_id : str
    job identifier from condor_submit stdout

  TODO
  ----
  are there python bindings for Condor that are more stable than parsing stdout?
  """
  pass

def get_revision_number():
  """
  Get revision number from git

  Raises
  ------
  CalledProcessError
    if git command fails to find the revision number this function will propagate the error:
    we do not handle this exception because reproducibility is too often overlooked

  ValueError
    if environment variable CS799_REPO_DIR is not set
  """
  env_var = 'CS799_REPO_DIR'
  repo_dir = os.environ.get(env_var)
  if repo_dir is None:
    raise ValueError("Missing environment variable {}".format(env_var))
  rev_number = sp.check_output(["get_revision_no.sh"])
  return rev_number.rstrip()

def log_script(argv):
  """
  Parameters
  ----------
  argv : list of str
    return value of sys.argv()
  """
  sys.stderr.write(" ".join(argv) + "\n")
  sys.stderr.write("Revision: {}\n".format(get_revision_number()))

def mkdir_p(dirpath):
  try:
    os.makedirs(dirpath)
  except OSError as e:
    pass
