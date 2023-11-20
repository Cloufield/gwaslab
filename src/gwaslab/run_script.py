import subprocess

def run_shell_script_single_command(script):
    plink_process = subprocess.Popen("exec script", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,text=True)
    output1,output2 = plink_process.communicate()
    plink_log+=output1 + output2+ "\n"
    plink_process.kill()

        try:
            output = subprocess.check_output(script, stderr=subprocess.STDOUT, shell=True, text=True)
            log.write(" -Saved results for CHR {} to : {}".format(i,"{}.clumps".format(out_single_chr)))
            plink_log +=output + "\n"
        except subprocess.CalledProcessError as e:
            log.write(e.output)