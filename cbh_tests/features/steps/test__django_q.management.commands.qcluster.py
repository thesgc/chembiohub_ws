from behave import given, when, then
import subprocess
import time
from signal import SIGINT
import select
import os

@given("I start the qcluster")
def step(context):
    """Use the exec function to ensure that the process can be killed with its subprocesses"""
    context.qfilename = "qclusterprocesslog.log"
    if os.path.exists(context.qfilename):
        os.remove(context.qfilename)
    cmd = "exec python manage.py qcluster"
    context.logfile = open(context.qfilename, "w")
    context.logfile.write("starting up\n")
    from django.conf import settings 

    context.django_q_process = subprocess.Popen(cmd, stderr=context.logfile, stdout=context.logfile, shell=True)
    counter = 0
    time.sleep(1)
    while True:
        context.logfile.flush()
        counter += 1
        if(counter > 150):
            raise Exception("Qcluster did not start")
        with open(context.qfilename, "r") as b:
            if  "running." in b.read():
                break
            else:
                print (b.read())

        time.sleep(1)





@when("I stop the qcluster")
def step(context):
    context.django_q_process.send_signal(SIGINT)
    context.django_q_process.wait()
    context.logfile.flush()
    time.sleep(1)
    context.logfile.close()
    


@then(u'I see the right qcluster output')
def step_impl(context):
    with open(context.qfilename, "r") as b:
        if  "has stopped." in b.read():
            pass
        else:
            raise Exception("qcluster did not stop")




