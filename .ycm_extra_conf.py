import os

def abspath(x):
  return os.path.abspath(os.path.normpath(os.path.expanduser(x)))

def copath(x):
  return os.path.abspath(os.path.normpath(os.path.join(__file__, "..", os.path.expanduser(x))))

base = os.path.normpath(os.path.realpath(os.path.dirname(__file__)))

flags = [
  '-Wall',
  '-Wextra',
  '-Werror',
  '-Wno-long-long',
  '-Wno-variadic-macros',
  '-fexceptions',
  '-std=c++14',
  '-x', 'c++',

  # ubuntu 16
  '-isystem', '/usr/include/c++/7',
  '-isystem', '/usr/include/x86_64-linux-gnu/c++/7',
  '-isystem', '/usr/include/c++/7/backward',
  '-isystem', '/usr/local/include',
  '-isystem', '/usr/include/x86_64-linux-gnu',
  '-isystem', '/usr/include',

  # freebsd (w/ clang-devel)
  '-isystem', '/usr/include/c++/v1',
  '-isystem', '/usr/local/llvm-devel/lib/clang/8.0.0/include',

  # mac
  '-isystem', '/opt/local/libexec/llvm-7.0/include/c++/v1',
  '-isystem', '/opt/local/libexec/llvm-7.0/lib/clang/7.0.1/include',
  '-isystem', '/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include',

  # zdock/libpdb
  '-I', 'contrib/eigen',
  '-I', 'contrib/Catch2/single_include',
  '-I', 'src/libpdb++',
  '-I', 'src/zdock',
  '-I', 'src/common',
  '-I', 'src/pdb',
  '-I', 'src',
  '-I', 'include',
]

def parsepaths(flags):
  stash = []
  for x in flags:
    if x not in ['-I', '-isystem']:
      if stash:
        if os.path.isdir(abspath(x)):
          while stash:
            yield stash.pop()
          yield abspath(x)
        elif os.path.isdir(copath(x)):
          while stash:
            yield stash.pop()
          yield copath(x)
        else:
          stash = []
      else:
        yield x
    else:
      stash.append(x)

def FlagsForFile(filename, *args, **kwargs):
  ret = {
    'flags': parsepaths(flags + ['-I', abspath(os.path.dirname(filename))]),
    'do_cache': True
  }
  return ret

