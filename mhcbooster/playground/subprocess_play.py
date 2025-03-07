import time
import multiprocessing
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map


def task(x):
    """Simulate a task and log progress in real-time."""
    time.sleep(2)  # Simulate processing time
    # print(f"Task {x} completed")  # Print logs directly
    return x * x

def worker(task_input):
    """Worker function to execute a task."""
    return task(task_input)

def run_parallel_tasks(tasks):
    """Run parallel tasks and use tqdm for progress bar."""
    # Use a Pool to execute tasks in parallel
    process_map(worker, tasks, max_workers=20, chunksize=1)
    # with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    #     # Use tqdm to display a progress bar for tasks
    #     results = list(tqdm(pool.imap(worker, tasks), total=len(tasks)))

# List of tasks to process
tasks = list(range(100))

# Run parallel tasks with progress
run_parallel_tasks(tasks)

