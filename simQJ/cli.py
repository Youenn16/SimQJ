import click
import subprocess
import os


params_ising_benchmark = [
    "ModelParameters:",
    "name=isingChain_density",
    "gamma=0.3",
    "h=0.7",
    "kappa=1.0",
    "J=1.0", 
    "",
    "SimulationParameters:",
    "L=8",
    "dt=0.005",
    "time_max=10",
    "no_click=False",
    "high_order=True",
    "record_probas=True",
    "",
    "DataMonitoring:",
    "time",
    "norm",
    "entanglement(2,4)",
    "density(0,1)",
]


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CPP_PATH = os.path.join(BASE_DIR, "bin", "engine_bin")

@click.group(chain=True)
def entry_point():
    """ SimQJ is a simulation toolbox for Stochastic Quantum Jumps trajectories ! """
    pass

@entry_point.command()
@click.option("--model", type=str, default="IsingModel", help="Model to benchmark")
def benchmark(model):
    """Benchmark in low dimension code with exact diagonalisation"""
    # click.secho("--- Starting C++ Module ---", fg="blue")

    with open(BASE_DIR+"/config/bench_config.param", "w") as f:
        for val in params_ising_benchmark:
            f.write(f"{val}\n")


    click.secho("--- Starting C++ Simulation ---", fg="blue")  
    result = subprocess.run( [BASE_DIR+"/build/simQJ", "-c", "config/bench_config.param"], capture_output=True, text=True )
    click.secho("--- Starting ED Python Simulation ---", fg="blue") 
    result = subprocess.run( ["python3", BASE_DIR+"/toolbox/exact_QJ.py", "--model", model], capture_output=True, text=True )
    click.secho("--- Plotting Comparison ---", fg="blue") 
    
# @click.option("--number", type=int, default=0, help="Specific size/site to plot")
# @click.argument("quantity", type=str, default="entanglement", required=True)

@entry_point.command(context_settings=dict(allow_interspersed_args=True))
@click.option("-q","--quantity", default="entanglement")
@click.option("-n","--number", type=int, default=-1)  
def plot(quantity, number):
    """Plot specific quantity of one trajectory""" 
    click.echo(f"Plotting {quantity} with number={number}")
    result = subprocess.run( ["python3", BASE_DIR+"/toolbox/visualize_traj.py", quantity, "-n", str(number)], capture_output=True, text=True ) 
    click.secho(result.stdout, fg="magenta") 
    click.secho(result.stderr, fg="magenta") 


 

# @entry_point.command()
# def command2():
#     """Runs the Python Script"""
#     click.secho("--- Starting Python Module ---", fg="green")
#     # You can import your other python functions here and call them directly
#     print("Python logic executed.")
    # click.secho(result.stdout, fg="magenta") 
    # click.secho(result.stderr, fg="red") 

if __name__ == "__main__":
    entry_point()