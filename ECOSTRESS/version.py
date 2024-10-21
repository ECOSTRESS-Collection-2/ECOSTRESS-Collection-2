from os.path import join, abspath, dirname


def version():
    with open(join(abspath(dirname(__file__)), 'version.txt')) as f:
        return f.read().strip()


def main():
    print(version())
