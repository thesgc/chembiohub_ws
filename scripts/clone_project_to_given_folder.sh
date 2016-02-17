set -e
FOLDER="$1"
sudo mkdir -p "$FOLDER"
sudo chown -R ubuntu "$FOLDER"
cd "$FOLDER"
sudo apt-get install -y git
git clone https://github.com/thesgc/chembiohub_ws .
git submodule init
git submodule update
git submodule foreach git pull origin master