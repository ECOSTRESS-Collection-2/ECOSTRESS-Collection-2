import subprocess
from os.path import abspath, dirname

directory = abspath(dirname(__file__))

def main():
    command = f'cd "{directory}" && julia --project=@. instantiate.jl'
    print(command)
    # system(command)
    subprocess.run(command, shell=True)

if __name__ == "__main":
    main()
