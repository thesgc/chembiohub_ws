
Feature: XLS upload
    I wish to submit an xls file that is either flagged as blinded / not blinded to a project I have editor rights. The XLS file either does /does not contain structural information. The substances in the XLS file are either pre-registered or not registered publically on the CBH resistration system.

	Scenario: User submits an XLS file containing unregistered substances to proja. The XLS file contains structural information and is flagged as not blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file contains structural data
        and I flag the file as not blinded
        and I submit file 
        then the XLS substances will be registered to proja

	Scenario: User submits an XLS file containing unregistered substances to proja. The XLS file contains structural information and is flagged as blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file contains structural data
        and I flag the file as blinded
        and I submit file 
        then the response would be Error - This XLS file contains structural information but has been flagged as blinded. Please flag as 'not blinded'.

    Scenario: User submits an XLS file containing unregistered substances to proja. The XLS file does not contain structural information and is flagged as blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file does not contain structural data
        and I flag the file as blinded
        and I submit file 
        then the XLS substances will be registered to proja


    Scenario: User submits an XLS file containing unregistered substances to proja. The XLS file does not contain structural information and is flagged as not blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file does not contain structural data
        and I flag the file as not blinded
        and I submit file 
        then the response would be Error - This XLS file contains no structural information but has been flagged as not blinded. Please flag as 'blinded'. 

Scenario: User submits an XLS file containing registered substances to proja. The XLS file contains structural information and is flagged as not blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file contains structural data
        and I flag the file as not blinded
        and I submit file 
        Then the response will be This substance has already been registered. Would you like to register as a new batch (or force registration)

	Scenario: User submits an XLS file containing registered substances to proja. The XLS file contains structural information and is flagged as blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file contains structural data
        and I flag the file as blinded
        and I submit file 
        then the responses would be This substance has already been registered. Would you like to register as a new batch (or force registration)? and Error - This XLS file contains structural information but has been flagged as blinded. Please flag as 'not blinded'. 

    Scenario: User submits an XLS file containing registered substances to proja. The XLS file does not contain structural information and is flagged as blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file does not contain structural data
        and I flag the file as blinded
        and I submit file 
        then the responses would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?>


    Scenario: User submits an XLS file containing registered substances to proja. The XLS file does not contain structural information and is flagged as not blinded.
	    Given a User
        When I log in 
        and a valid project exists proja
        and I have editor rights for proja
        when I submit an XLS file
        and a substance in the XLS file is not pre-registered
        and the file does not contain structural data
        and I flag the file as not blinded
        and I submit file 
        then the response would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?> and <Error - This XLS file does not contain structural information but has been flagged as not blinded. Please flag as 'blinded'.>








        