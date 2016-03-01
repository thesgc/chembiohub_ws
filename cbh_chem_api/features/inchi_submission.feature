Feature: INCHI upload

'''
INCHI submission

Feature: INCHI upload
    I wish to submit substances to a project I have editor rights using multiple INCHIs. The substances are either pre-registered or not registered publically on the CBH resistration system. The INCHIs entered are either all dissimilar or contain duplications


	Scenario: User submits unregistered substances to proja in a INCHI format. The INCHIs are all dissimilar to one another. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's INCHIs are not pre-registered
        and the INCHIs are disimilar
        and I have editor rights for proja
        when I submit INCHIs
        then the substances will be registered to proja


    Scenario: User submits unregistered substances to proja in a INCHI format. Some of the INCHIs are duplicated. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's INCHIs are not pre-registered
        and some of the INCHIs are duplicated
        and I have editor rights for proja
        when I submit INCHIs
        then the the response would be <Some of your INCHIs are duplicated. Please remove duplications and resubmit> 


    Scenario: User submits registered substances to proja in a INCHI format. Some of the INCHIs are duplicated. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's INCHIs are pre-registered
        and some of the INCHIs are duplicated
        and I have editor rights for proja
        when I submit INCHIs
        then the the response would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?> and <Some of your INCHIs are duplicated. Please remove duplications and resubmit>   


    Scenario: User submits registered substances to proja in a INCHI format. The INCHIs are all dissimilar to one another.
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's INCHIs are pre-registered
        and the INCHIs are all dissimilar.
        and I have editor rights for proja
        when I submit INCHIs
        then the the response would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?>  

'''