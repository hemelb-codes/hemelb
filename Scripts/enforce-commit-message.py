import re

# Use via
# [hooks]
# pretxncommit = python:/path_to_hemelb/Scripts/enforce-commit-message.py:checkCommitMessage

def checkCommitMessage(ui, repo, **kwargs):
    """
    Checks a single commit message for adherence to commit message rules.  

    To use add the following to your project .hg/hgrc for each
    project you want to check, or to your user hgrc to apply to all projects.

    [hooks]
    pretxncommit = python:path/to/script/enforce-message.py:checkCommitMessage
    """

    hg_commit_message = repo['tip'].description()

    if(checkMessage(hg_commit_message) == True):
        ui.warn('This is not a valid commit message: ' + hg_commit_message + '\n')
        printUsage(ui)
        return True
    else:
        return False
        
def checkMessage(message):
    if not re.match('\A\W*merge\W*\Z', message, re.IGNORECASE) == None:
        return False
    if not re.match('\A.*\#\d+.*\Z', message) == None:
        return False
    return True

def printUsage(ui):
    ui.warn('Commit message must be the single word \'Merge\' or contain \'#<ticket_id>\'\n')
