#!/usr/bin/env python3
# Copyright (C) 2021 FIUBioRG
# SPDX-License-Identifier: MIT
#
# -*-python-*-

import logging
import os
from os import environ, getenv
from os.path import abspath, realpath, dirname
import sys
import subprocess

if sys.version_info.major < 3:
    logging.error("Python3 is required to run this SConstruct")

print("Updating csv-parser...")
subprocess.call(['git', 'submodule', 'update', '--init'])

AddOption(
    "--pluma",
    dest="pluma",
    type="string",
    nargs=1,
    action="store",
    metavar="DIR",
    default=abspath("../../src"),
    help="DIR for the installation of PluMA's header files"
)

env = Environment(
    Env=environ,
    CC=getenv("CC", "cc"),
    CXX=getenv("CXX", "c++"),
    SHCCFLAGS=["-fpermissive", "-fPIC", "-I.", "-O2"],
    SHCXXFLAGS=["-std=c++11", "-fpermissive", "-fPIC", "-I.", "-O2"],
    CCFLAGS=["-fpermissive", "-fPIC", "-I.", "-O2"],
    CXXFLAGS=["-std=c++11", "-fpermissive", "-fPIC", "-I.", "-O2"],
    LICENSE=["MIT"],
    SHLIBPREFIX="lib"
)

env.Append(
    SHCXXFLAGS=["-I" + env.GetOption("pluma")],
    CXXFLAGS=["-I" + env.GetOption("pluma")]
)

if not sys.platform.startswith("darwin"):
    env.Append(LINKFLAGS=["-rdynamic"])
    env.Append(LIBS=["rt"])
else:
    env.Append(CCFLAGS=["-DAPPLE"])

config = Configure(env)

if not config.CheckCC():
    Exit(1)

if not config.CheckCXX():
    Exit(1)

if not config.CheckSHCXX():
    Exit(1)

if not config.CheckCXXHeader("csv_parser/csv_parser.hpp"):
    Exit(1)

libs = [
    "m",
    "c"
]

for lib in libs:
    if not config.CheckLib(lib):
        Exit(1)

plugin = realpath("./ATriaPlugin.cpp")
filesInPath = Glob(dirname(plugin) + "./*.cpp")
pluginName = plugin.replace(".cpp", "")

env.SharedLibrary(
    target=pluginName,
    source=filesInPath
)