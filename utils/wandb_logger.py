import wandb

class wandb_project():
    def __init__(self, project_name, config):
        self.project_name = project_name
        self.config = config
        
        self.run = wandb.init(
            project=self.project_name,
            config=config
        )

    def log(self, metrics):
        self.run.log(
            metrics
        )


