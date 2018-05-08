
import utils.file_utils as utils
import fungcat_settings as f_settings


def get_uniprot_to_species():
    uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1, 2)
    return uniprot_to_species


def setupEmailCommand(recipient, subject=None, message=None, from_wyatt=False):
    """ Build a command for sending an email using an SMTP server with the 'sendEmail' package
    Currently very specific to my (Jeff's) setup
    *from_wyatt* if this is not being sent from wyatt (False), then ssh to wyatt to send the email
    """
    # default message if there were no errors
    if message is None:
        message = "Finished `date`"
    #print "Sending email to %s" % (recipient)
    # currently sent from the vm I have setup on wyatt
    # using my personal gmail (with a special app password)
    # TODO also send an email if there was an error
    email_command = "sendemail -f jeffreynlaw@gmail.com -t %s " % (recipient) + \
            "-u \"%s\" -m \"%s\" " % (subject, message) + \
            "-s smtp.gmail.com:587 -o tls=yes -xu jeffreynlaw@gmail.com -xp \"lxkiidlmmjpekfie\""
    # if this command is going to be run on a different computer, then ssh to wyatt first
    if from_wyatt is False:
        email_command = "ssh -t jeffl@wyatt.cs.vt.edu \"%s\"" % (email_command)
    return email_command
