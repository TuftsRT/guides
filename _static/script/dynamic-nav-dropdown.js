// Dynamically move navbar items in/out of dropdown based on available space
(function () {
  function getDropdownBtn() {
    const navbar = document.querySelector(".bd-navbar-elements");
    if (!navbar) return;
    const dropdownBtn = navbar.querySelector(".dropdown");
    if (!dropdownBtn) {
      const btn = document.createElement("li");
      btn.className = "nav-item dropdown";
      btn.innerHTML = `
        <button class="btn dropdown-toggle nav-item" type="button"
          data-bs-toggle="dropdown" aria-expanded="false"
          aria-controls="pst-nav-more-links">
        More
        </button>
        <ul id="pst-nav-more-links" class="dropdown-menu">
        </ul>
      `;
      btn.style.display = "none";
      if (document.documentElement.dataset.theme === "dark") {
        btn.querySelector(".dropdown-menu").classList.add("dropdown-menu-dark");
      }
      navbar.appendChild(btn);
      return btn;
    }
    return dropdownBtn;
  }

  function getCopyWidth(element) {
    if (!element) return 0;
    const clone = element.cloneNode(true);
    clone.style.visibility = "hidden";
    clone.style.position = "absolute";
    clone.style.display = "block";
    clone.style.left = "-9999px";
    clone.style.top = "0";
    document.body.appendChild(clone);
    const width = clone.offsetWidth;
    document.body.removeChild(clone);
    return width;
  }

  function getNavbarWidth(navbarContainer, navItems) {
    navItems.forEach((item) => (item.style.display = "none"));
    let navbarWidth = navbarContainer.offsetWidth;
    navItems.forEach((item) => {
      item.removeAttribute("style");
    });
    return navbarWidth;
  }

  function getHorizontalMargins(element) {
    if (!element) return 0;
    const style = window.getComputedStyle(element);
    const left = parseFloat(style.marginLeft) || 0;
    const right = parseFloat(style.marginRight) || 0;
    return left + right;
  }

  function getElements() {
    const navbarContainer = document.querySelector(
      ".navbar-header-items__center"
    );
    const navbar =
      navbarContainer && navbarContainer.querySelector(".bd-navbar-elements");
    const dropdownMenu = navbar && navbar.querySelector(".dropdown-menu");
    const dropdownItems =
      dropdownMenu && Array.from(dropdownMenu.querySelectorAll(":scope > li"));
    const navItems =
      navbar &&
      Array.from(navbar.querySelectorAll(":scope > li:not(.dropdown)"));
    return { navbarContainer, navbar, dropdownMenu, navItems, dropdownItems };
  }

  function countItemsThatFit(elements, margin, availableWidth) {
    let count = 0;
    let usedWidth = 0;
    for (let i = 0; i < elements.length; i++) {
      let itemWidth = getCopyWidth(elements[i]) + margin;
      if (usedWidth + itemWidth > availableWidth) break;
      usedWidth += itemWidth;
      count++;
    }
    return count;
  }

  function moveDropdownItemToNavbar(item, navbar, dropdownBtn) {
    if (!item || !navbar || !dropdownBtn) return;
    item.classList.add("nav-item");
    const link = item.querySelector("a");
    if (link) {
      link.classList.remove("dropdown-item");
    }
    navbar.insertBefore(item, dropdownBtn);
  }

  function moveNavbarItemToDropdown(item, dropdownMenu) {
    if (!item || !dropdownMenu) return;
    item.classList.remove("nav-item");
    const link = item.querySelector("a");
    if (link) {
      link.classList.add("dropdown-item");
    }
    dropdownMenu.insertBefore(item, dropdownMenu.firstChild);
  }

  function updateNavbar(dropdownBtn, dropdownBtnWidth, margin) {
    const { navbarContainer, navbar, dropdownMenu, navItems, dropdownItems } =
      getElements();
    if (
      !navbarContainer ||
      !navbar ||
      !dropdownBtn ||
      !dropdownMenu ||
      !navItems ||
      !dropdownItems
    )
      return;
    let navbarWidth = getNavbarWidth(navbarContainer, navItems);
    let navItemsWidth =
      navItems.reduce((total, item) => total + item.offsetWidth, 0) +
      margin * navItems.length;
    let dropdownItemsWidth =
      dropdownItems.reduce((total, item) => total + getCopyWidth(item), 0) +
      margin * dropdownItems.length;
    if (navbarWidth >= navItemsWidth + dropdownItemsWidth) {
      dropdownItems.forEach((item) => {
        moveDropdownItemToNavbar(item, navbar, dropdownBtn);
      });
      dropdownBtn.style.display = "none";
      document
        .querySelector(".navbar-header-items__center")
        .removeAttribute("style");
    } else if (navbarWidth >= navItemsWidth + dropdownBtnWidth) {
      let originalNumberOfItems = navItems.length;
      let toMove = countItemsThatFit(
        dropdownItems,
        margin,
        navbarWidth - navItemsWidth - dropdownBtnWidth
      );
      for (let i = 0; i < toMove; i++) {
        moveDropdownItemToNavbar(dropdownItems.shift(), navbar, dropdownBtn);
      }
      if (originalNumberOfItems === 0 && toMove > 0) {
        dropdownBtn.querySelector("button").textContent = " More ";
        document
          .querySelector(".navbar-header-items__center")
          .removeAttribute("style");
      }
    } else if (navbarWidth < navItemsWidth + dropdownBtnWidth) {
      let toMove =
        navItems.length -
        countItemsThatFit(navItems, margin, navbarWidth - dropdownBtnWidth);
      for (let i = 0; i < toMove; i++) {
        moveNavbarItemToDropdown(navItems.pop(), dropdownMenu);
      }

      dropdownBtn.removeAttribute("style");
      if (navItems.length === 0 && toMove > 0) {
        dropdownBtn.querySelector("button").textContent = " Menu ";
        document.querySelector(
          ".navbar-header-items__center"
        ).style.justifyContent = "right";
      }
    }
  }

  let margin = 0;
  let dropdownBtn;
  let dropdownBtnWidth;
  window.addEventListener("DOMContentLoaded", () => {
    dropdownBtn = getDropdownBtn();
    dropdownBtnWidth = getCopyWidth(dropdownBtn);
    margin = getHorizontalMargins(dropdownBtn);
    updateNavbar(dropdownBtn, dropdownBtnWidth, margin);
  });
  window.addEventListener("resize", () =>
    updateNavbar(dropdownBtn, dropdownBtnWidth, margin)
  );
})();
