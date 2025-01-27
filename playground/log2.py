import sys
import time


# Custom class to redirect stdout to both the console and a file with a timestamp
class DualOutput:
    def __init__(self, file):
        self.file = file
        self.console = sys.stdout  # Save original stdout (console)

    def write(self, message):
        if message != '\n':
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            message = f"{timestamp} - {message}"

        # Write the message to both console and file
        self.console.write(message)
        self.file.write(message)

    def flush(self):
        # Flush the output (this is required for compatibility)
        self.console.flush()
        self.file.flush()


# A sample function that generates output
def my_function():
    print("This is a standard output message.")  # Prints to console and log file
    time.sleep(50)
    raise Exception("This is an error message.")  # Will raise an exception


# Function to run the Python function and capture output to a log file
def log_function_output():
    # Open the log file in write mode
    with open('output.log', 'w') as log_file:
        # Redirect sys.stdout to both console and log file
        sys.stdout = DualOutput(log_file)

        try:
            # Run the function
            my_function()

        except Exception as e:
            # Log the exception to both console and log file
            print(f"Caught an exception: {e}")

        finally:
            # Restore the original stdout
            sys.stdout = sys.__stdout__


# Run the function and log its output
log_function_output()
