(function () {
  document.addEventListener("DOMContentLoaded", function () {
  // Bot icon SVG — speech bubble with a simple robot face
  var icon = `<svg xmlns="http://www.w3.org/2000/svg" width="30" height="30" viewBox="0 0 64 64" aria-hidden="true" focusable="false">
    <!-- rounded head -->
    <rect x="10" y="14" width="44" height="32" rx="8" fill="#fff"/>
    <!-- antenna -->
    <line x1="32" y1="14" x2="32" y2="6" stroke="#fff" stroke-width="3" stroke-linecap="round"/>
    <circle cx="32" cy="5" r="3" fill="#fff"/>
    <!-- eyes -->
    <rect x="18" y="25" width="10" height="8" rx="3" fill="#4AABE3"/>
    <rect x="36" y="25" width="10" height="8" rx="3" fill="#4AABE3"/>
    <!-- mouth -->
    <rect x="20" y="37" width="24" height="4" rx="2" fill="#4AABE3"/>
    <!-- chin triangle / speech pointer -->
    <polygon points="24,46 40,46 32,54" fill="#fff"/>
  </svg>`;

  var link = document.createElement("a");
  link.id = "chatbot-bubble";
  link.href = "/hpc/chat.html";
  link.setAttribute("aria-label", "Ask our Tufts Research Technology Guides AI Assistant");

  var btn = document.createElement("span");
  btn.id = "chatbot-bubble-btn";
  btn.innerHTML = icon;

  var label = document.createElement("span");
  label.id = "chatbot-bubble-label";
  label.textContent = "Ask our Tufts Research Technology Guides AI Assistant";

  // Label appears to the left of the button
  link.appendChild(label);
  link.appendChild(btn);

  document.body.appendChild(link);
  });
})();
