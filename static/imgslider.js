let picSlideIndex = 1;
showPicSlides(picSlideIndex);

function picPlusSlides(n) {
  showPicSlides(picSlideIndex += n);
}

function picCurrentSlide(n) {
  showPicSlides(picSlideIndex = n);
}

function showPicSlides(n) {
  let i;
  let slides = document.getElementsByClassName("imgmySlides");
  let dots = document.getElementsByClassName("imgdot");
  
  if (n > slides.length) { picSlideIndex = 1; }
  if (n < 1) { picSlideIndex = slides.length; }
  
  for (i = 0; i < slides.length; i++) {
    slides[i].style.display = "none";
  }
  
  for (i = 0; i < dots.length; i++) {
    dots[i].className = dots[i].className.replace(" active", "");
  }
  
  slides[picSlideIndex-1].style.display = "block";
  dots[picSlideIndex-1].className += " active";
}

setInterval(() => {
  picPlusSlides(1);
}, 5000);