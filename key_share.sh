#!/bin/bash
# https://linuxcluster.wordpress.com/2015/04/27/simple-bash-script-to-setup-shared-ssh-keys-on-cluster/

# Exit script on Error
set -e

# Check for SSH Directory
if [ ! -d ~/.ssh ]; then
   mkdir -p ~/.ssh/
fi


# Check for existence of passphrase
if [ ! -f ~/.ssh/id_rsa.pub ]; then
        ssh-keygen -t rsa -N "" -f ~/.ssh/id_rsa
        echo "Execute ssh-keygen --[done]"
fi

# Check for existence of authorized_keys and append the shared ssh keys
if [ ! -f ~/.ssh/authorized_keys ]; then
        touch ~/.ssh/authorized_keys
        echo "Create ~/.ssh/authorized_keys --[done]"
        chmod 700 ~/.ssh/authorized_keys
        cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
        echo "Append the public keys id_rsa into authorized keys --[done]"
        chmod 400 ~/.ssh/authorized_keys
        chmod 700 ~/.ssh/
fi

# Create user's ssh config it not exist
if [ ! -f ~/.ssh/config ]; then
        touch ~/.ssh/config
        echo "StrictHostKeyChecking no" > ~/.ssh/config
        echo "StrictHostKeyChecking no --[done]"
        chmod 644 ~/.ssh/config
fi
# Unset error on exit or it will affect after bash command 🙂
set +e

