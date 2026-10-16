#!/bin/env bash
set -e

# bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_links.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_links.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_signalp.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_deep.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_effectorp.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_antismash.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_iprscan.sh
bash /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_annotate.sh
# /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_compare.sh
# /home/aubgxs001/comparative_fungi/scripts/funannotate_links/runMe_funannotate_orthofinder.sh