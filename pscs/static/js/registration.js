// This file is intended to be used with the register.html template, which is
// used for the creation of new accounts.

function isPasswordMatch(){
  let password = document.getElementById('password');
  let confirmation = document.getElementById('passwordConfirm');
  return password.value === confirmation.value;
}

function validatePassword(){
  if( isPasswordMatch()){
    return true;
  }
  else{
    let confirmation = document.getElementById("passwordConfirm");
    confirmation.setCustomValidity("Password mismatch");
    confirmation.reportValidity();
    confirmation.setCustomValidity("");
    return false;
  }
}

function validateUsername(){
  let isValid = true;
  let usernameEl = document.getElementById('username');
  let username = usernameEl.value;
  let firstChar = username.substr(0,1);
  if(username.length === 0){  // check for empty field
    usernameEl.setCustomValidity('Username is required.');
    isValid = false;
  }
  else if(firstChar.toLowerCase() === firstChar.toUpperCase()){  // returns true iff first char is not alphabetical
    usernameEl.setCustomValidity('First character in username must be alphabetical');
    isValid = false;
  }
  if(!isValid){
    usernameEl.reportValidity();
    usernameEl.setCustomValidity("");
  }
  return isValid;
}

function validatePHI(){
  let checkEl = document.getElementById('noPHI');
  if(!checkEl.checked){
    checkEl.setCustomValidity("Agreement is necessary for registration.");
    checkEl.reportValidity();
    checkEl.setCustomValidity("");
  }

  return checkEl.checked;
}

function validateDataUse(){
  let checkEl = document.getElementById("dataUse");
  if(!checkEl.checked){
    checkEl.setCustomValidity("Agreement is necessary for registration.");
    checkEl.reportValidity();
    checkEl.setCustomValidity("");
  }
  return checkEl.checked;
}

function validateForm(){
  let isValid = true;
  // user validate
  isValid &&= validateUsername();
  // Password
  isValid &&= validatePassword();
  // email validated by input; need to pass to server to check against db
  // PHI ack
  isValid &&= validatePHI();
  // data use ack
  isValid &&= validateDataUse();

  return isValid;
}