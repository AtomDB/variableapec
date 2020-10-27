from setuptools import setup
# from https://gist.github.com/dcreager/300803 with "-dirty" support added
import os, Extension

linapprox =  Extension("linear_approx",['linear_approx.c'],\
                       define_macros = [('MAJOR_VERSION', '1'),\
                                        ('MINOR_VERSION','0')])

def GetGitVersion():
    '''report the git commit/branch/tag on which we are '''
    repo = Repo(".", search_parent_directories=True)
    git = repo.git    

    branchOrTag=git.rev_parse('--abbrev-ref', 'HEAD')

    if branchOrTag == 'HEAD':
        # get tag name
        # since several tags are sometime put on the same commit we want to retrieve all of them
        # and use the last one as reference
        # Note:
        # branchOrTag=`git describe --tags --exact-match` does not provided the latest created tag in case several point to the same place
        currentSha=git.rev_parse('--verify','HEAD')

        # list all tags on the current sha using most recent first:
        allTags=git.tag('--points-at',currentSha,'--sort=-creatordate')
        print (allTags)

        allTagsArray=allTags.split(' ') #create an array (assuming tags are separated by space)

        # if we checkouted a commit with no tag associated, the allTagsArray is empty we can use directly the sha value
        if len(allTagsArray) == 0:
            branchOrTag=git.rev-rev_parse('--short','HEAD') # take the short sha
        else:

            branchOrTag=allTagsArray[0] #first from the list
    else:
        #add the head commit id on the current branch
        branchOrTag="{}[{}]".format(branchOrTag,git.rev_parse('--short', 'HEAD'))

    return branchOrTag

if on_rtd:
  from mock import Mock as MagicMock

  extmos= []

  class Mock(MagicMock):
      @classmethod
      def __getattr__(cls, name):
              return Mock()

  MOCK_MODULES = ['liblinapprox']
  sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

else:
  extmos= [linapprox]

setup(name='variableapec',
      version=variableapec.__version__,
      description='variableapec tool for AtomDB python library.',
      url='http://www.atomdb.org',
      author='Keri Heuer',
      author_email='kh3286@drexel.edu',
      license='Smithsonian',
      classifiers=['Development Status :: 4 - Beta',\
                   'Environment :: Console',\
                   'Intended Audience :: Developers',\
                   'Intended Audience :: Education',\
                   'Intended Audience :: End Users/Desktop',\
                   'Intended Audience :: Science/Research',\
                   'Topic :: Scientific/Engineering :: Astronomy',\
                   'Topic :: Scientific/Engineering :: Physics',\
                   'Programming Language :: Python :: 3',\
                   'Operating System :: POSIX'],
      zip_safe=False,
      long_description = README_TEXT,\
      install_requires=[
      "requests",\
      "wget",\
      "numpy>=1.9.0",\
      "scipy>=1.4.0",\
      "joblib",\
      "mock",\
      "astropy",\
      "pycurl"],
      ext_modules = extmos)
