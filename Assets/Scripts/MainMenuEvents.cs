using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;
public class MainMenuEvents : MonoBehaviour
{
    private UIDocument _document;

    private Button _button;

	private List<Button> _menubuttons = new List<Button>();
    private void Start()
    {
        _document = GetComponent<UIDocument>();

       // _button = _document.rootVisualElement.Q("EnterButton") as Button;
	   // _button.RegisterCallback<ClickEvent>(OnPlayGameClick);

		_menubuttons = _document.rootVisualElement.Query<Button>().ToList();

		for (int i = 0; i < _menubuttons.Count; i++)
		{
			_menubuttons[i].RegisterCallback<ClickEvent>(AllButtonsClick);
		}
	}

	

	private void OnDisable()
	{
        //_button.UnregisterCallback<ClickEvent>(OnPlayGameClick);

		for (int i = 0; i < _menubuttons.Count; i++)
		{
			_menubuttons[i].UnregisterCallback<ClickEvent>(AllButtonsClick);
		}
	}

	/*private void OnPlayGameClick(ClickEvent evt)
    {
        Debug.Log("Kill Yourself");
    }*/

	private void AllButtonsClick(ClickEvent evt)
	{
		Debug.Log("Test");
	}
}