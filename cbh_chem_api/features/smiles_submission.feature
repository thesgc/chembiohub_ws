
Feature: SMILES upload
    I wish to submit substances to a project I have editor rights using multiple SMILES. The substances are either pre-registered or not registered publically on the CBH resistration system. The SMILES entered are either all dissimilar or contain duplications


	Scenario: User submits unregistered substances to proja in a SMILES format. The SMILES are all dissimilar to one another. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's SMILES are not pre-registered
        and the SMILES are disimilar
        and I have editor rights for proja
        when I submit SMILES
        then the substances will be registered to proja


    Scenario: User submits unregistered substances to proja in a SMILES format. Some of the SMILES are duplicated. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's SMILES are not pre-registered
        and some of the SMILES are duplicated
        and I have editor rights for proja
        when I submit SMILES
        then the the response would be <Some of your SMILES are duplicated. Please remove duplications and resubmit> 


    Scenario: User submits registered substances to proja in a SMILES format. Some of the SMILES are duplicated. 
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's SMILES are pre-registered
        and some of the SMILES are duplicated
        and I have editor rights for proja
        when I submit SMILES
        then the the response would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?> and <Some of your SMILES are duplicated. Please remove duplications and resubmit>   


    Scenario: User submits registered substances to proja in a SMILES format. The SMILES are all dissimilar to one another.
	    Given a User
        When I log in 
        and a valid project exists proja
        and the submitted substance's SMILES are pre-registered
        and the SMILES are all dissimilar.
        and I have editor rights for proja
        when I submit SMILES
        then the the response would be <This substance has already been registered. Would you like to register as a new batch (or force registration)?>    