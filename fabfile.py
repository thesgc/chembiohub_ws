from fabric.operations import prompt, local, run, sudo
from fabric.api import cd, prefix

from fabric.state import env

from fabric.context_managers import shell_env



try:
    from server_configs import SERVER_CONFIGS
    
except ImportError:
    print """ You Must provide a file server_configs.py to use fabric for deployment.
    It should have a format:

SERVER_CONFIGS = [{"host_string": "my-server.com", 
                    "user" : "astretton",
                    "directory": "/srv/chembiohub/chembiohub_ws", 
                    "conda_env_path": "/home/chembiohub/anaconda2/envs/hub",
                    "url": "https://my-server.com/hub",
                    "user_running_chembiohub" : "chembiohub"},
                    {"host_string": "my-server2.com", 
                    "user" : "astretton",
                    "directory": "/srv/chembiohub/chembiohub_ws", 
                    "conda_env_path": "/home/chembiohub/anaconda2/envs/hub",
                    "url": "https://my-server2.com/hub",
                    "user_running_chembiohub" : "chembiohub"}]
"""


def deploy():
    for index, conf in enumerate(SERVER_CONFIGS):
        print ("%d: %s" % (index + 1, conf["url"]) ) 
        
    prompt("Which instance of the application would you like to deploy to?", key="instance_number")
    instance = int(env.instance_number) -1
    env.update(SERVER_CONFIGS[instance])
    new_path = env.conda_env_path + "/bin:$PATH"
    with shell_env(PATH=new_path, CONDA_ENV_PATH=env.conda_env_path, ):
        with cd(env.directory):
            sudo("git pull origin master", user=env.user_running_chembiohub)
            sudo("git submodule foreach git pull origin master", user=env.user_running_chembiohub)
           # sudo("conda install -y --file anaconda_requirements.txt" , user=env.user_running_chembiohub)
            sudo("pip install -r pip_requirements.txt", user=env.user_running_chembiohub)
            sudo("python manage.py migrate", user=env.user_running_chembiohub)

        with cd(env.directory + "/src/ng-chem"):
            sudo("bower install", user=env.user_running_chembiohub)
        with cd(env.directory):
            sudo("python manage.py collectstatic", user=env.user_running_chembiohub)
        sudo("supervisorctl reload")
