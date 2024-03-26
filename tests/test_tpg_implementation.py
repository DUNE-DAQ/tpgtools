import pytest
import os
import subprocess
import time

from pathlib import Path
from tempfile import TemporaryDirectory

patterns = ["patt_golden", "patt_square", "patt_square_left", "patt_square_right"]
#patterns = ["patt_golden", "patt_square"]

channels = [
  pytest.param(x) for x in [str(_) for _ in range(0, 64, 1)]
  ]

clock_ticks = [
  pytest.param(x) for x in [str(_) for _ in range(1, 64, 1)]
  ]


@pytest.fixture(scope="session")
def tmp_dir():
  with TemporaryDirectory() as tmp_dir:
    tmp_path = Path(tmp_dir)

@pytest.fixture(scope="session")
def tmp_dir_path(tmp_path_factory):
  path = tmp_path_factory.mktemp("tpgtest")
  yield path


"""
# -------------------------------------------------------------
# TESTS
"""

@pytest.fixture(params=patterns)
def get_pattern(request):
  yield request.param

@pytest.fixture(params=channels)
def get_channel(request):
  yield request.param

@pytest.fixture(params=clock_ticks)
def get_tick(request):
  yield request.param

@pytest.fixture(params=["10","64","499"])
def get_threshold(request):
  yield request.param

@pytest.fixture(params=["1","2"])
def get_frames(request):
  yield request.param

@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--printout")


@pytest.fixture()
def get_input(tmp_dir_path,get_pattern,get_channel,get_tick,get_threshold,get_frames, cmdopt):
  if cmdopt == "more":
    print("\nDBG start : ", ) 
    print("DBG start path : ", tmp_dir_path)

  tick = get_tick
  pattern = get_pattern
  channel = get_channel
  threshold = get_threshold
  num_frames = get_frames
  print("\nDBG start cmdopt : ", cmdopt)
  if cmdopt == "more":
    print("DBG start tick: ", tick)
    print("DBG start pattern: ", pattern)
    print("DBG start channel: ", channel)
    print("DBG start threshold: ", threshold)
    print("DBG start number of frames: ", num_frames)

  str_path = os.path.join(*[tmp_dir_path,pattern,"chan_"+channel,"tick_"+tick])
  path = Path(str_path) 
  subprocess.call(["mkdir -p " + str_path], shell=True)

  # pattern generation
  exe = "wibeth_tpg_pattern_generator "
  args =  ["-f", "/cvmfs/dunedaq.opensciencegrid.org/assets/files/d/d/1/wibeth_output_all_zeros.bin"]
  args += ["-o", str_path]
  args += ["--save-trigprim -w"]
  args += ["-n", num_frames]
  args += ["-t", threshold]
  args += ["-i", channel]
  args += ["-c", tick]
  args += ["-p", pattern]
  args += ["-s", "__"+tick]
  args = ' '.join(args)
  cmd =  exe + args
  if cmdopt == "more":
    print("\nDBG command: ", cmd)

  output_frame_file = os.path.join(*[str_path,pattern+"_chan_"+channel+"_tick_"+tick+"_wibeth_output.bin"])
  if not os.path.exists(output_frame_file):
    fp = path / str("log_pattgen__" + tick + ".txt");
    f = open(str(fp), "w")
    subprocess.call([cmd], shell=True, stdout=f)
  if cmdopt == "more":
    subprocess.call(["ls -ltr " + str_path], shell=True)

  yield (pattern, channel, tick, threshold, num_frames, path, str_path, output_frame_file, cmdopt)

@pytest.mark.parametrize('alg', ["SimpleThreshold", "AbsRS"])
def test_all_params(get_input, alg):

  (pattern, channel, tick, threshold, num_frames, path, str_path, output_frame_file, cmdopt) = get_input

  str_subpath = os.path.join(*[str_path,alg])
  subpath = Path(str_subpath)
  subprocess.call(["mkdir -p " + str_subpath], shell=True)

  if cmdopt == "more":
    print("\nDBG finish pattern: ", pattern)
    print("DBG finish channel: ", channel)
    print("DBG finish tick: ", tick)
    print("DBG finish threshold: ", threshold)
    print("DBG finish num_frames: ", num_frames)
    print("DBG finish path: ", path)
    print("DBG finish str_path: ", str_path)
    print("DBG finish algorithm: ", alg)

  # validation
  fp1 = subpath / str("log_naive__" + tick + ".txt");
  fp2 = subpath / str("log_avx__"+ tick + ".txt");
  f1 = open(str(fp1), "w")
  f2 = open(str(fp2), "w")
  exe = " ".join(("wibeth_tpg_workload_emulator", ""))
  args =  ["-f", output_frame_file]
  args += ["-o", str_subpath]
  args += ["-r false -m VDColdboxChannelMap --save-trigprim"]
  args += ["-n", num_frames]
  args += ["-t", threshold]
  args += ["-c", tick]
  args += ["-a", alg]
  args = ' '.join(args)
  cmd =  exe + args
  impl_naive = " ".join(("-i", "NAIVE", ""))
  impl_avx = " ".join(("-i", "AVX", ""))

  cmd1 = cmd + " " + impl_naive
  cmd2 = cmd + " " + impl_avx
  if cmdopt == "more":
    print("\nDBG cmd1: ", cmd1)
    print("\nDBG cmd2: ", cmd2)

  subprocess.call([cmd1], shell=True, stdout=f1)
  subprocess.call([cmd2], shell=True, stdout=f2)
  if cmdopt == "more":
    print("\nDBG str_subpath ", str_subpath)
    subprocess.call(["ls -ltr " + str_path], shell=True)
    subprocess.call(["ls -ltr " + str_subpath], shell=True)
  subprocess.call(["mv " + str_subpath + "/TP_dump_AVX*"   + " " + str_subpath + "/naive.txt"], shell=True)
  subprocess.call(["mv " + str_subpath + "/TP_dump_NAIVE*" + " " + str_subpath + "/avx.txt"], shell=True)
  if cmdopt == "more":
    subprocess.call(["ls -ltr " + str_path], shell=True)
    subprocess.call(["ls -ltr " + str_subpath], shell=True)

  fp1 = subpath / "naive.txt";
  fp2 = subpath / "avx.txt";
  assert fp1.read_text() == fp2.read_text()

# ------------------------------------------------------------------------

