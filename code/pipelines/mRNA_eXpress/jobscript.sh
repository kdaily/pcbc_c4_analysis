#!/bin/sh
#rule: {job}
#input: {job.input}
#output: {job.output}
cd $SGE_O_WORKDIR
source /home/ubuntu/.bashrc
workon snakemake
snakemake --snakefile {self.workflow.snakefile} \
--force -j{self.cores} \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
> /dev/null && touch "{jobfinished}" || touch "{jobfailed}"
exit 0
