// This file is intended to be used with the register.html template, which is
// used for the creation of new accounts.

function isPasswordMatch(){
  password = document.getElementById('password');
  confirmation = document.getElementById('passwordConfirm');
  console.log(password.value);
  console.log(confirmation.value);
  if(password.value != confirmation.value){
    return false;
  }
  return true;
}

function validatePassword(){
  if( isPasswordMatch()){
    return true;
  }
  else{
    passwordConfirm.setCustomValidity("Password mismatch");
    passwordConfirm.reportValidity();
    return false;
  }
}

function validateUsername(){
  isValid = true;
  usernameEl = document.getElementById('username');
  username = usernameEl.value;
  firstChar = username.substr(0,1);
  console.log(firstChar);
  console.log(username);
//  alert(firstChar);
  if(username.length == 0){  // check for empty field
    usernameEl.setCustomValidity('Username is required.');
    isValid = false;
  }
  else if(firstChar.toLowerCase() == firstChar.toUpperCase()){  // returns true iff first char is not alphabetical
    usernameEl.setCustomValidity('First character in username must be alphabetical');
    isValid = false;
  }
  usernameEl.reportValidity();
  return isValid;
}

function validatePHI(){
  checkEl = document.getElementById('noPHI');
  return checkEl.checked;
}

function validateForm(){
  isValid = true;
  // user validate
  isValid &&= validateUsername();
  // Password
  isValid &&= validatePassword();
  // email
  // already validated by input
  // PHI ack
  isValid &&= validatePHI();
  return isValid;
}