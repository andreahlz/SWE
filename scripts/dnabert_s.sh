#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pr_root=${script_dir%/scripts}

source "$pr_root/config.env"
 
output_dir=mkdir -p "$pr_root/$out_dir_dnaberts"


