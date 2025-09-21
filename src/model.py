from runner.runner import Runner

if __name__ == "__main__":

    run = (
        Runner(run_type = "model")
            .parse_arguments()
            .create_run_directories()
            .copy_config_data()
            .submit_abaqus_command()
    )
