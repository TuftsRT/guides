(function () {
  document.addEventListener("DOMContentLoaded", function () {
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

    var message = "Ask our Tufts Research Technology Guides AI Assistant";
    
    var wrapper = document.createElement("div");
    wrapper.id = "chatbot-bubble-wrap";

    var link = document.createElement("a");
    link.id = "chatbot-bubble";
    link.href = "/hpc/chat.html";
    link.setAttribute("aria-label", message);

    var btn = document.createElement("span");
    btn.id = "chatbot-bubble-btn";
    btn.innerHTML = icon;

    var label = document.createElement("span");
    label.id = "chatbot-bubble-label";
    label.textContent = message;

    link.appendChild(label);
    link.appendChild(btn);

    // Mobile popup — shown on first tap instead of navigating directly
    var popup = document.createElement("div");
    popup.id = "chatbot-bubble-popup";
    popup.setAttribute("hidden", "");

    var popupText = document.createElement("p");
    popupText.textContent = message;

    var popupLink = document.createElement("a");
    popupLink.id = "chatbot-popup-link";
    popupLink.href = "/hpc/chat.html";
    popupLink.textContent = "Open Chat";

    popup.appendChild(popupText);
    popup.appendChild(popupLink);

    wrapper.appendChild(popup);
    wrapper.appendChild(link);

    var isTouchDevice = function () {
      return window.matchMedia("(hover: none) and (pointer: coarse)").matches;
    };

    link.addEventListener("click", function (e) {
      if (isTouchDevice()) {
        e.preventDefault();
        if (popup.hasAttribute("hidden")) {
          popup.removeAttribute("hidden");
          link.setAttribute("aria-expanded", "true");
          popupLink.focus();
        } else {
          popup.setAttribute("hidden", "");
          link.setAttribute("aria-expanded", "false");
        }
      }
    });

    document.addEventListener("click", function (e) {
      if (!wrapper.contains(e.target)) {
        popup.setAttribute("hidden", "");
        link.setAttribute("aria-expanded", "false");
      }
    });

    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && !popup.hasAttribute("hidden")) {
        popup.setAttribute("hidden", "");
        link.setAttribute("aria-expanded", "false");
        link.focus();
      }
    });

    var footer = document.querySelector("footer.bd-footer");
    if (footer) {
      document.body.insertBefore(wrapper, footer);
    } else {
      document.body.appendChild(wrapper);
    }
  });
})();
