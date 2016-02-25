set -e
FOLDER="$1"
USERNAME="$2"
echo "$FOLDER"
sudo mkdir -p "$FOLDER"
sudo chown -R "$USERNAME" "$FOLDER"
cd "$FOLDER"
sudo apt-get install -y git
git clone https://github.com/thesgc/chembiohub_ws .
git submodule init
git submodule update
git submodule foreach git pull origin master