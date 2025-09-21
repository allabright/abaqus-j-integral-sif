from runner.runner import Runner

if __name__ == "__main__":

    run = (
        Runner(run_type = "export")
            .parse_arguments()
            .create_run_directories()
            .copy_and_read_input_data()
            .submit_abaqus_command()
    )
